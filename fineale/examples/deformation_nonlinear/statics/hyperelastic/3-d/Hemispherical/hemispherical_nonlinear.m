% Nonlinear twisted beam
%
function [lambdas,u1s,u2s]=hemispherical_nonlinear
E=6.825e7;
nu=0.3;
thickness = 0.04;
ncirc=1; nlayers=2;
% analytical solution for the vertical deflection under the load
analyt_sol=0.0924;
graphics = true;
R=10;
utol = thickness/10000;
graphics = true;
nincr  =10;
tolerance = thickness/100;

% prop = property_deformation_neohookean (struct('E',E,'nu',nu));
% mater = material_deformation_neohookean_triax(struct('property',prop));
prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_stvk_triax(struct('property',prop));

surface_integration_rule=gauss_rule(struct('dim',2,'order',2));


%          Mesh
%% Create the mesh and initialize the geometry
[fens,fes]=H_spherical_shell_mesh(R-thickness/2,thickness,ncirc,nlayers,[]);
%  [fens,fes]=H_spherical_shell_18deg_hole_mesh(R-thickness/2,thickness,ncirc,ncirc,nlayers, []);

femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
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
component= [2];
fixed_value= 0;
node_list = fenode_select (fens,struct ('box',[-10000 10000 0 0 -10000 10000],'inflate',tolerance));
u= set_ebc(u,node_list,true,component,fixed_value);
component= [1];
fixed_value= 0;
node_list = fenode_select (fens,struct ('box',[0 0 -10000 10000 -10000 10000],'inflate',tolerance));
u= set_ebc(u,node_list,true,component,fixed_value);
component= [];
fixed_value= 0;
node_list = fenode_select (fens,struct ('box',[0 0 0 0 R-thickness/2 R-thickness/2],'inflate',tolerance));
u= set_ebc(u,node_list,true,component,fixed_value);

u   = apply_ebc (u);

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 100;
u = u*0; % zero out the displacement
utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= [-0.171332705794298  -7.882139089645855   5.594516542362686   4.394378011771107  -1.931989037461593   1.264389523440495                   0   0   1.000000000000000  54.988185976473318];
end

node_list1 =fenode_select (fens,struct ('box',[-inf inf 0 0 0  0],'inflate',tolerance));
force1= [1/length(node_list1);0;0];
node_list2 =fenode_select (fens,struct ('box',[0 0 -inf inf 0  0],'inflate',tolerance));
force2= [0;-1/length(node_list2);0];
nffemm1 = femm_deformation_linear (struct ('material',[],...
    'fes',fe_set_P1(struct('conn',reshape(node_list1,[],1))),...
    'integration_rule',point_rule));
nffemm2 = femm_deformation_linear (struct ('material',[],...
    'fes',fe_set_P1(struct('conn',reshape(node_list2,[],1))),...
    'integration_rule',point_rule));

u1s=[0];u2s=[0];lambdas=[0];
femm  =associate_geometry(femm,geom); 
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
        
        fi= force_intensity(struct('magn',t*force1));
        FL = distrib_loads(nffemm1, sysvec_assembler, geom, 0*u, fi, 0);
        fi= force_intensity(struct('magn',t*force2));
        FL = FL + distrib_loads(nffemm2, sysvec_assembler, geom, 0*u, fi, 0);
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
    [~,femm] = restoring_force(femm,sysvec_assembler,geom,u1,u);        % final update
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
    us{end+1} =u;
    u1=mean(gather_values(u,node_list1));
    u2=mean(gather_values(u,node_list2));
    u1s=[u1s,mean(u1(:,1))];u2s=[u2s,mean(u2(:,2))];lambdas=[lambdas,t];
    incr = incr + 1;
end

end
