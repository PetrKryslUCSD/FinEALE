% Free-free beam with self equilibrated loads
%
function free_free_beam
E=1000;
nu=0.249;
rho= 1.0 ;
W=2.5;
H=5;
L= 25;
htol=min([L,H,W])/1000;
magn = -0.2*12.2334/4;
Force =magn*W*H*2;
Force*L^3/(3*E*W*H^3*2/12);
graphics = ~true;
stress_range = [-30,30];
u_scale=1;
nW=2;nL=5;nH=2;

% Material
prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
mater = material_deformation_linear_triax (struct('property',prop ));
%     Mesh it
% [fens,fes]= T10_block(W,L,H,nW,nL,nH);
% [fens,fes]= T4_blocka(W,L,H,nW,nL,nH);
[fens,fes]= H8_block(W,L,H,nW,nL,nH);
fens=translate_mesh(fens,-[W,L,H]/2);
                
%     Make the finite element model machine
% femm = femm_deformation_linear (struct ('material',mater, 'fes',fes,...
%     'integration_rule',tet_rule (struct('npts',5))));
femm = femm_deformation_linear (struct ('material',mater, 'fes',fes,...
    'integration_rule',gauss_rule (struct('dim',3,'order',2))));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out

% Number equations
u   = numberdofs (u);
% Assemble the system matrix
K = stiffness(femm, sysmat_assembler_sparse, geom, u);

bdry_fes = mesh_boundary(fes, []);
bclp = fe_select(fens, bdry_fes, ...
    struct ('box',[-inf inf L/2 L/2 -inf inf],'inflate',htol));
bfesp=femm_deformation_linear (struct ('material',mater, 'fes',subset(bdry_fes,bclp),...
    'integration_rule',gauss_rule (struct('dim',2,'order',2))));
bclm = fe_select(fens, bdry_fes, ...
    struct ('box',[-inf inf -L/2 -L/2 -inf inf],'inflate',htol));
bfesm=femm_deformation_linear (struct ('material',mater, 'fes',subset(bdry_fes,bclm),...
    'integration_rule',gauss_rule (struct('dim',2,'order',2))));

F = zeros(u.nfreedofs,1);

% Twisting torque
fi= force_intensity(struct('magn',@(x)[+x(3);0;-x(1)]));
F = F + distrib_loads(bfesp,sysvec_assembler, geom, u, fi, 2);
fi= force_intensity(struct('magn',@(x)-[+x(3);0;-x(1)]));
F = F + distrib_loads(bfesm,sysvec_assembler, geom, u, fi, 2);
% Bending torque
% fi= force_intensity(struct('magn',@(x)[0;+x(3);0]));
% F = F + distrib_loads(bfesp,sysvec_assembler, geom, u, fi, 2);
% fi= force_intensity(struct('magn',@(x)-[0;+x(3);0]));
% F = F + distrib_loads(bfesm,sysvec_assembler, geom, u, fi, 2);

[V,D]=eig(full(K));
[~,ix]= sort(diag(real(D)));
D=  real(D(ix,ix));
V=real(V(:,ix));

for j=1:6
K = K + 1e-6*V(:,j)*V(:,j)';
end 
% Solve
u = scatter_sysvec(u, K\F);

% Plot
gv=graphic_viewer;
gv=reset (gv,struct ([]));
scale=10;
cmap = jet;
draw(femm,gv, struct ('x', geom, 'u', +scale*u,'shrink',1.0));
draw(femm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1));

set_graphics_defaults(gcf)

end