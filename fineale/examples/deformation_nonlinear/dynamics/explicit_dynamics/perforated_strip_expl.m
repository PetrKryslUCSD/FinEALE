function perforated_strip_expl
graphics= true;
 
L= 10; % Length of the plate
W = 1; % Thickness of the plate
H = 10; % Width
R= 3;
nL=4;nH=4;nR=3;
magn=60;
tapp=0.05;
tend =  2*tapp;
rho=1e-3;
lambda = 0.75e4;
G= 50;

nu =1/2*lambda/(G+lambda)
E=G*(3*lambda+2*G)/(G+lambda);


[fens,fes]=Q4_elliphole(R,R,L,H,nL,nH,nR,[]);
[fens,fes] = H8_extrude_Q4(fens,fes,3,@(x,i)([x,0]+[0,0,W/3*i]));
% [fens,fes] = H8_to_H27(fens,fes);

prop = property_deformation_linear_iso (struct('E',E,'nu',nu,'rho',rho));
mater = material_deformation_stvk_triax(struct('property',prop));
femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero outfunction five
% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_comp=ones(1,length (ebc_fenids))*1;
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,0,0,-Inf,Inf],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_comp=ones(1,length (ebc_fenids))*2;
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_comp=ones(1,length (ebc_fenids))*3;
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);
% loaded boundary
bdry_fes = mesh_boundary(fes, []);
bcl = fe_select(fens, bdry_fes, ...
    struct ('box',[L,L,-Inf,Inf,-Inf,Inf],'inflate',H/1000));
efemm = femm_deformation_nonlinear(struct ('material',mater, 'fes',subset(bdry_fes,bcl),...
    'integration_rule',gauss_rule(struct('dim',2,'order',4))));

M =lumped_mass(femm, sysmat_assembler_sparse, geom, u);
K = stiffness(femm, sysmat_assembler_sparse, geom, u,u);

    o2=eigs(K,M,1,'LM');
    dt= 0.99* 2/sqrt(o2) / 20 % compensate for serious reduction in thickness


if graphics,
gv=graphic_viewer;
    gv=clear(gv,[]);
    gv=reset(gv,[]);
end


%     Velocity field
v=0*u;

% Solve
[femm] = update(femm,  geom, u, u);
t=0;
U0= gather_sysvec(u);
V0= gather_sysvec(v);
fi= force_intensity(struct ('magn',[magn;0; 0;]));
FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
F1 = FL + restoring_force(femm,sysvec_assembler, geom,u,u);       % Internal forces
A1=M\(F1);
A0 =A1;
step =0;
while t <tend
    step = step  +1;
     t=t+dt;
    fi= force_intensity(struct ('magn',((t<tapp)*(t/tapp)+(t>=tapp)*1.0)*[magn;0; 0;]));
    FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
    % Update displacement
    U1 = U0 +dt*V0+(dt^2)/2*A0;% displacement update
    u = scatter_sysvec(u, U1);
    F1 = FL + restoring_force(femm,sysvec_assembler, geom,u,u);       % Internal forces
    % Compute the new acceleration.
    A1=M\(F1);
    % Update the velocity
    V1 = V0 +(dt/2)* (A0+A1);
    % Bring the the displacement and velocity fields up to date
    v = scatter_sysvec(v, V1);
    if graphics && (mod(step,100) == 0)
        disp(['Time ' num2str(t)])
        gv=reset(gv,[]);    draw(femm,gv, struct ('x', geom,'u',u, 'facecolor','none'));
        figure (gcf); pause(0.05); cam =camget(gv);
    end
    % Switch the temporary vectors for the next step.
    U0 = U1;
    V0 = V1;
    A0 = A1;
    if (t==tend)
        break;
    end
    if (t+dt>tend)
        dt =tend-t;
    end
    t=t+dt;
end



end


