% 45/-45 two-layer plate
%
function [u1s,lambdas]=two_layer
pu= physical_units_struct;
W=100*pu.MM;
L=400*pu.MM;
ts= [5,5]*pu.MM;
nl=16; nt=2*[1,1]; nw=6;
nl=32; nt=4*[1,1]; nw=12;
nl=8; nt=1*[1,1]; nw=3;
p=  10e3*pu.PA;
utol = min(ts)/1e10;
gtol = min(ts)/1e5;
graphics = true;
nincr  =10;

% Anisotropic 
E1=5*1e6*pu.PA; E2=1*1e6*pu.PA; G12=0.5e6*pu.PA; nu= 0.25;
angles =[+45,-45];
    function Rm = LayerRm(XYZ, ts, label)% label equals the layer number here
        Rm= rotmat(angles(label)/180*pi* [0,0,1]);
    end

prop = property_deformation_linear_transv_iso(...
    struct('E1',E1,'E2',E2,'G12',G12,'nu12',nu,'nu23',nu));
mater = material_deformation_bb_transviso_triax(struct('property',prop));


surface_integration_rule=gauss_rule(struct('dim',2,'order',4));


%          Mesh
%% Create the mesh and initialize the geometry
[fens,fes] = H8_composite_plate(L,W,ts,nl,nw,nt);;
[fens,fes] = H8_to_H64(fens,fes);


femm = femm_deformation_nonlinear(struct (...
    'material',mater,'Rm',@LayerRm, ...
    'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',4))));

% Material
% Finite element block
bdry_fes = mesh_boundary(fes, []);
bcl = fe_select(fens, bdry_fes, ...
    struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',gtol));

efemm = femm_deformation (struct ('material',mater, 'fes',subset(bdry_fes,bcl),...
    'integration_rule',surface_integration_rule));
sfemm = femm_deformation (struct ('material',mater, 'fes',bdry_fes,...
    'integration_rule',surface_integration_rule));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
fixed_value= 0;
node_list = fenode_select (fens,struct ('box',[0 0 -100*W 100*W -100*W 100*W],'inflate',gtol));;
u= set_ebc(u,node_list,true,1,fixed_value);
u= set_ebc(u,node_list,true,2,fixed_value);
u= set_ebc(u,node_list,true,3,fixed_value);
u   = apply_ebc (u);

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 10;
u = u*0; % zero out the displacement
utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= [-0.171332705794298  -7.882139089645855   5.594516542362686   4.394378011771107  -1.931989037461593   1.264389523440495                   0   0   1.000000000000000  54.988185976473318];
end

enl=fenode_select (fens,struct ('box',[L L 0 0 0 0],'inflate',gtol));


u1s=[]; lambdas=[];
incr=1;
while (incr <= nincr)
    t = incr* tup / nincr;
    disp(['Increment ' num2str(incr) ]); % pause
    % Initialization
    u1 = u; % guess
    u1 = apply_ebc(u1);
    du = 0*u; % this will hold displacement increment
    du = apply_ebc(du);
    iter=1;
    %             gv=reset(clear(gv,[]),[]);
    %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
    femm1=femm;
    while 1
        Load=zeros(3,1); Load(1) =p*t;
        fi=force_intensity(struct('magn',Load));
        FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
        %         sum(FL)
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
        K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
        % Displacement increment
        du = scatter_sysvec(du, K\F);
        R0 = dot(F,gather_sysvec(du));
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1+du,u);       % Internal forces
        R1 = dot(F,gather_sysvec(du));
        a = R0/R1;
        if ( a<0 )
            eta = a/2 +sqrt((a/2)^2 -a);
        else
            eta =a/2;
        end
        if (imag(eta)~=0)
            disp('######################  Inverted elements?')
        end
        eta=min( [eta, 1.0] );
        u1 = u1 + eta*du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(abs(du.values)) < utol) break; end;                    % convergence check
        iter=iter+1;
    end
    [ignore,femm] = restoring_force(femm,sysvec_assembler,geom,u1,u);        % final update
    disp(['    Converged for t=' num2str(t)]); % pause
    u = u1;                                               % update the displacement
    if graphics
        gv=reset(clear(gv,[]),[]);
        draw(sfemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
        draw(sfemm,gv, struct ('x', geom, 'u', u,'facecolor','y', 'shrink',1.0));
        %         camset (gv,Cam);
        interact(gv);
        pause(0.5); %Cam =camget(gv);
    end
    us{end+1} =u;
    uend=gather_values(u,enl);
    u1s=[u1s,uend];
    lambdas=[lambdas,t];
    incr = incr + 1;
end
u1s  =reshape(u1s',[3,nincr])';
end
