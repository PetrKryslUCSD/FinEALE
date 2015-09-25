% Square domain with internal heat source  and fixed temperature 
% on the boundary.
function pnpSquare
%% 
% Evaluate for a given data
a=2.5; dy=a/2*sin(15/180*pi); dx=a/2*cos(15/180*pi); Q=4.5; k=1.8; Dz=1.0;
x=[0,0; 2*dx,-2*dy; a,0; 2*dx,2*dy];  % Coordinates
[fens,fes]=Q4_quadrilateral(x,1,1,struct('other_dimension' ,1))
prop=property_heat_diffusion(struct('thermal_conductivity',k,'source',Q));
mater=material_heat_diffusion (struct('property',prop));
femm = femm_heat_diffusion (struct ('material',mater,...
    'fes',fes,...
    'integration_rule',gauss_rule(struct('dim', 2, 'order',2))));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta = nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
theta = set_ebc(theta, [2,3,4], true, 1, 0.0);
theta = apply_ebc (theta);
theta = numberdofs (theta);
K = conductivity(femm, sysmat_assembler_sparse, geom, theta);
F = zeros(theta.nfreedofs,1);
fi= force_intensity(struct('magn',Q));
F = F + distrib_loads(femm, sysvec_assembler, geom, theta, fi, 3);
F = F + nz_ebc_loads_conductivity(femm, sysvec_assembler, geom, theta);
theta = scatter_sysvec(theta, K\F);
% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',hot));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, T)));
geomT=nodal_field(struct ('name', ['geomT'], ...
    'data',[geom.values, theta.values]));
draw(fes, gv, struct ('x',geomT, 'u',0*geomT,...
        'colorfield',colorfield, 'shrink',0.9));
    draw(fes, gv, struct ('x',geom, 'u',0*geom, ...
        'facecolor','none'));
     draw_integration_points(femm, gv, struct('x',geom,'theta',theta, 'scale', 0.1))

axis equal vis3d
set_graphics_defaults(gcf)
labels