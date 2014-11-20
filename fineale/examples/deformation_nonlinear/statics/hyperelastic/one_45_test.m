function Perforated_strip_neohook_linesearch
u= physical_units_struct;
graphics= ~~true;
stabfact=0.1;
L= 100*u.MM; % Length of the plate
H = 50*u.MM; % Width
W = 10*u.MM; % Thickness of the plate
nL=1;nH=1;nW=1;
magn=9e6*u.PA;
maxit = 15;
nincr =10;
utol = 0.000001;

    % Anisotropic, but much stronger anisotropy
    E1=1000*1e6*u.PA; E2=1000*1e6*u.PA; E3=1*1e6*u.PA; G12=0.2e6*u.PA;  G13=G12; G23=0.2e6*u.PA;
    nu12= 0.25; nu13= 0.25; nu23= 0.25;
    aangle =+45;
Rm= rotmat(aangle/180*pi* [0,0,1]);
    % [xestim, beta, c, residual] = richextrapol([],[1,2,4])
    prop = property_deformation_linear_ortho (...
        struct('E1',E1,'E2',E2,'E3',E3,...
        'G12',G12,'G13',G13,'G23',G23,...
        'nu12',nu12,'nu13',nu13,'nu23',nu23));
    stabprop = property_deformation_linear_iso (struct('E',E3,'nu',0));
    

[fens,fes] = H8_block(L,H,W,nL,nH,nW);
% [fens,fes] = H8_to_H27(fens,fes);

mater = material_deformation_stvk_triax(struct('property',prop));
stabmater = material_deformation_stvk_triax(struct('property',stabprop));
femm = femm_deformation_nonlinear_h8msgs(struct ('material',mater,...
    'stabilization_material',stabmater, 'Rm',Rm,'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero outfunction five
% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, 1, ebc_val);
ebc_fenids=fenode_select (fens,struct('box',[0,0,0,0,-Inf,Inf],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, 2, ebc_val);
ebc_fenids=fenode_select (fens,struct('box',[0,0,0,0,0,0],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, 3, ebc_val);
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
        if (max(abs(du.values)) < utol) break; end;% convergence check
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
    draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','none','alpha', 0.2));
end
%         camset(gv, [  -75.6602  613.9688  210.2074   20.6651   38.3517 20.5774         0         0    1.0000    5.3084]);
clear K
end

