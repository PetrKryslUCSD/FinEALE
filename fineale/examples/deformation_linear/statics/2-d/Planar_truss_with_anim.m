% Example 18: planar truss, example 3.7, p.72 from Hutton
disp('Example: planar truss, with animated deflection');

% Parameters:
E=1e7;

% Mesh: finite element nodes
fens= fenode_set(struct ('xyz',[[0 0 0]; ...
    [0 40 0];...
    [40 0 0];...
    [40 40 0];...
    [80 0 0];...
    [80 40 0]]));
%Mesh: finite elements
fes = fe_set_L2(struct (...
    'conn',[[1 3];[1 4];[2 4];[3 4];[3 5];[5 4];[6 4];[5 6]],...
    'other_dimension',  1.5));

ebc_fenids=[1 1 1 2 2 2];
ebc_prescribed=[1 1 1 1 1 1];
ebc_comp=[1 2 3 1 2 3];
ebc_val=ebc_comp*0;

% Material
prop = property_deformation_linear_iso(struct('E',E));
mater = material_deformation_linear_uniax (struct('property',prop));
% Finite element model machine
femm = femm_deformation_linear (struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('order',1,'dim', 1 )),...
    'Rm','geniso'));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
for i=1:count(fens) % this loop sets the zero out-of-plane displacement
    u=set_ebc(u, [i],[1],[3],[0]);
end
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);
% Assemble the stiffness matrix
K = stiffness(femm, sysmat_assembler_sparse,    geom, u);

% Load
n=force_intensity(struct('magn',[0;-2000;0]));
lfemm = femm_deformation_linear (struct ('fes',fe_set_P1(struct('conn',3)), ...
    'integration_rule',point_rule));
F = distrib_loads(lfemm, sysvec_assembler, geom, u, n, 0);
n=force_intensity(struct('magn',[+2000; 0;0]));
lfemm = femm_deformation_linear (struct ('fes',fe_set_P1(struct('conn',5)), ...
    'integration_rule',point_rule));
F = F + distrib_loads(lfemm, sysvec_assembler, geom, u, n, 0);
n=force_intensity(struct('magn',[+4000; +6000;0]));
lfemm = femm_deformation_linear (struct ('fes',fe_set_P1(struct('conn',6)), ...
    'integration_rule',point_rule));
F = F + distrib_loads(lfemm, sysvec_assembler, geom, u, n, 0);

% Solve
u = scatter_sysvec(u, K\F);
% get(u,'values')

% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
scale=100;
draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','blue'));
draw(femm,gv, struct ('x', geom, 'u', (scale)*u,'facecolor','red'));
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on

draw_axes(gv, struct('length',10 ));
view(2);

pause (2);
for scale=sin(-pi/2+(0:1:84)/21*2*pi)
    gv=reset (gv,struct ('limits',[-5 100 -5 80 -20 20]));
    draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','blue'));
    draw(femm,gv, struct ('x', geom, 'u', 100*(scale+1)*u,'facecolor','red'));
    pause(0.1);
end