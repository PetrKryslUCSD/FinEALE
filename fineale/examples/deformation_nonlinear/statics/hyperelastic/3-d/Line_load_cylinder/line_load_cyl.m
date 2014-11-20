
% function [u1s,lambdas]=line_load_cyl
function res=line_load_cyl
% Reese (2002)  On the equivalence of mixed....
mu=6000;
lam=24000;
% [x,y]=solve ('24000-x*y/( 1+y)/(1-2*y)','6000-x/2/(1+y)')
E=16800;
nu=2/5;
L = 15;
R = 9;
thickness = 2.0;
n=16; nt=1;
graphics = ~true;
utol = thickness/1e9;
nincr  =1;
gtol = thickness/100;

prop = property_deformation_neohookean (struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));
alp=min([1.0,max([(1-2*nu)*3,0.001])]);;
nuhat =1/2-(1-2*nu)/2/alp;%A reasonable choice.   No locking anyway.
stabprop = property_deformation_neohookean (struct('E',E,'nu',nuhat));
stabmater = material_deformation_neohookean_triax(struct('property',stabprop));

surface_integration_rule=gauss_rule(struct('dim',2,'order',2));


%          Mesh
%% Create the mesh and initialize the geometry
[fens,fes]= H8_block(180/360*2*pi,L,thickness,n,n/2,nt);
xyz=fens.xyz;
for i=1:count (fens)
    a=xyz(i,1); y=xyz(i,2); z=xyz(i,3);
    xyz(i,:)=[(R-thickness/2+z)*sin(a) y (R-thickness/2+z)*cos(a)];
end
fens.xyz=xyz;
% drawmesh({fens,fes},'fes');
            
femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater,...
    'fes',fes, 'stabilization_material',stabmater,...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));

% Material
% Finite element block
bdry_fes = mesh_boundary(fes, []);

sfemm = femm_deformation (struct ('material',mater, 'fes',bdry_fes,...
    'integration_rule',surface_integration_rule));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
component= [1];
node_list = fenode_select (fens,struct ('box',[0 0 -inf inf -inf inf],'inflate',gtol));
u= set_ebc(u,node_list,true,component,0.0);
% u   = apply_ebc (u);
% show_field_as_markers(graphic_viewer,struct('x',geom,'u',u,'nl',node_list))
component= [2];
node_list = fenode_select (fens,struct ('box',[-inf inf L L -inf inf],'inflate',gtol));
u= set_ebc(u,node_list,true,component,0.0);
% u   = apply_ebc (u);
% show_field_as_markers(graphic_viewer,struct('x',geom,'u',u,'nl',node_list))
component= [3];
node_list = fenode_select (fens,struct ('box',[0 0 -inf inf -R-thickness/2 -R-thickness/2],'inflate',gtol));
u= set_ebc(u,node_list,true,component,0.0);

u   = apply_ebc (u);

% show_field_as_markers(graphic_viewer,struct('x',geom,'u',u,'nl',node_list))

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 7500;
u = u*0; % zero out the displacement
utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= [-0.171332705794298  -7.882139089645855   5.594516542362686   4.394378011771107  -1.931989037461593   1.264389523440495                   0   0   1.000000000000000  54.988185976473318];
end

cl=fe_select(fens,bdry_fes,struct('box',[0 0 -inf inf 0 +R+thickness/2],'inflate',gtol));
ffes=subset(bdry_fes,cl);
bdry_ffes=mesh_boundary(ffes,[]);
cl=fe_select(fens,bdry_ffes,struct('box',[0 0 -inf inf +R+thickness/2 +R+thickness/2],'inflate',gtol));
linefes=subset(bdry_ffes,cl);

nffemm2 = femm_deformation_linear (struct ('material',[],...
    'fes',linefes,'integration_rule',trapezoidal_rule(struct('dim',1))));

lA = fenode_select (fens,struct ('box',[0 0 0 0 +R+thickness/2 +R+thickness/2],'inflate',gtol));

u1s=[0]; lambdas=[0];
femm  =update(femm,geom,u,u); 
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
        
        fi= force_intensity(struct('magn',t*[0;0;-1/15]));
        FL = distrib_loads(nffemm2, sysvec_assembler, geom, 0*u, fi, 1);
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
        %                 draw(sfemm,gv, struct ('x', geom, 'u', u1,'facecolor','none'));
        %                 %                 draw(sfemm,gv, struct ('x', geom, 'u', 0*u1,'facecolor','y'));
        %                 interact(gv); pause(1);
        if (max(abs(du.values)) < utol) break; end;                    % convergence check
        iter=iter+1;
    end
    femm  =update(femm,geom,u1,u);
    disp(['    Converged for t=' num2str(t)]); % pause
    u = u1;                                               % update the displacement
    if graphics
        gv=reset(clear(gv,[]),[]);
        draw(sfemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
        draw(sfemm,gv, struct ('x', geom, 'u', u,'facecolor','y', 'shrink',1.0));
        camset (gv,Cam);
        interact(gv);
        pause(0.5); Cam =camget(gv);
    end
    u1A=gather_values(u,lA);
    u1s=[u1s,u1A(3)]; lambdas=[lambdas,t];
    incr = incr + 1;
end

res=[n,u1s(end)];
end
