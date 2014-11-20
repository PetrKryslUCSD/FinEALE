function Perforated_strip_neohook_linesearch
graphics= ~~true;
stabfact=0.1;
L= 10; % Length of the plate
W = 1; % Thickness of the plate
H = 10; % Width
R= 3;
nL=4;nH=2;nR=2;
magn=90;
maxit = 15;
nincr =10;
utol = 0.0001;

lambda = 0.75e4;
G= 50;

nu =1/2*lambda/(G+lambda)
E=G*(3*lambda+2*G)/(G+lambda);


[fens,fes]=Q4_elliphole(R,R,L,H,nL,nH,nR,[]);
[fens,fes] = H8_extrude_Q4(fens,fes,1,@(x,i)([x,0]+[0,0,W*i]));
[fens,fes] = H8_to_H27(fens,fes);

prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_stvk_triax(struct('property',prop));
femm = femm_deformation_nonlinear(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',3))));
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

gv=graphic_viewer;
dt= 1/nincr;
t=0; % time=load magnitude
incr=1;
while (incr <= nincr)
    t=t+ dt;
    disp(['Increment ' num2str(incr) ]); % pause
    % Initialization
    u1 = u; % guess
    du = 0*u; % this will hold displacement increment
    u1 = apply_ebc(u1);
    du = apply_ebc(du);
    
    iter=1;
    if graphics,
        gv=clear(gv,[]);
        gv=reset(gv,[]);
        %             axis off;
        %             draw(femm,gv, struct ('x', geom,'u',u, 'facecolor','cyan'));
        %             draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','Yellow','alpha', 0.2));
        % %             camset(gv,cam);
        %             cam =camget(gv);
        %                  print(gcf, '-r50', [mfilename n2s0p(Frame) '.png'], '-dpng');
    end
    while 1
        femm1=femm;
        fi= force_intensity(struct ('magn',[t*magn;0; 0;]));
        FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
        K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
        du = scatter_sysvec(du, K\F); % Displacement increment
        R0 = dot(F,gather_sysvec(du));
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1+du,u);       % Internal forces
        R1 = dot(F,gather_sysvec(du));
        a = R0/R1;
        if ( a<0 )
            eta = a/2 +sqrt((a/2)^2 -a);
        else
            eta =a/2;
        end
        eta=min( [eta, 1.0] );
        u1 = u1 + eta*du;   % increment displacement
        if graphics,
            draw(femm1,gv, struct ('x', geom,'u',u1, 'facecolor','none'));
            figure (gcf); pause(0.5); cam =camget(gv);
        end
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(abs(eta*du.values)) < utol) break; end;% convergence check
        if (iter >maxit)% bailout for failed convergence
            error( [' Possible failed convergence']);
        end
        iter=iter+1;
    end
    [ignore,femm] = restoring_force(femm,sysvec_assembler,geom,u1,u);        % final update
    disp(['    Converged for t=' num2str(t)]); % pause
    u = u1;       % update the displacement
    incr = incr + 1;
end
if graphics,
    gv=reset(clear(gv,[]),[]);camset (gv,cam);
    draw(femm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
    draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','Yellow','alpha', 0.2));
end
%         camset(gv, [  -75.6602  613.9688  210.2074   20.6651   38.3517 20.5774         0         0    1.0000    5.3084]);
clear K
end

