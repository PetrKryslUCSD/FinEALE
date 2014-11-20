% Rectangular domain, heat conduction through the thickness of a wall.
function pnpWall4T3
pum=physical_units_machine;
kappa1=0.1*pum('W/K/m'); % thermal conductivity, material 1
kappa2=0.3*pum('W/K/m'); % thermal conductivity, material 2
Ta1=(273.15+10)*pum('K'); % ambient temperature on the left-hand side surface of the wall
h1=5*pum('W/K/m^2');% surface heat transfer coefficient
Ts2=(273.15-20)*pum('K'); % temperature of the right-hand side surface of the wall
L=0.15*pum('m');% thickness of the wall
fens=fenode_set();
fens.xyz=[0,0; L/2,0; L,0; 0,L/3; L/2,L/3; L,L/3];
fes1T3 = fe_set_T3(struct('conn',[4,1,5; 2,5,1], 'other_dimension',1.0));
fes2T3 = fe_set_T3(struct('conn',[5,2,6; 3,6,2], 'other_dimension',1.0));
fes1L2 = fe_set_L2(struct('conn',[4,1], 'other_dimension',1.0));

mater1=material_heat_diffusion (struct(...
    'property',property_heat_diffusion(struct('thermal_conductivity',kappa1*eye(2)))));
mater2=material_heat_diffusion (struct(...
    'property',property_heat_diffusion(struct('thermal_conductivity',kappa2*eye(2)))));
femm1T3 = femm_heat_diffusion (struct ('material',mater1,...
    'fes',fes1T3,...
    'integration_rule',tri_rule(struct('npts', 1))));
femm2T3 = femm_heat_diffusion (struct ('material',mater2,...
    'fes',fes2T3,...
    'integration_rule',tri_rule(struct('npts', 1))));
femm1L2 = femm_heat_diffusion (struct ('material',mater1,...
    'fes',fes1L2,'surface_transfer',h1,...
    'integration_rule',gauss_rule(struct('dim',1,'order', 2))));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta = nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
fenids=[fenode_select(fens,struct('box',[L L 0 L/3],...
    'inflate', 0.01))];
fixed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+Ts2;
theta = set_ebc(theta, fenids, fixed, comp, val);
theta = apply_ebc (theta);
theta = numberdofs (theta);
K = conductivity(femm1T3, sysmat_assembler_sparse, geom, theta)...
    + conductivity(femm2T3, sysmat_assembler_sparse, geom, theta);
H= surface_transfer(femm1L2, sysmat_assembler_sparse, geom, theta);
fenids= connected_nodes(femm1L2.fes);
amb = set_ebc(0*theta, fenids, true, [], Ta1);
amb = apply_ebc (amb);
F = nz_ebc_loads_conductivity(femm2T3, sysvec_assembler, geom, theta)...
    + surface_transfer_loads(femm1L2, sysvec_assembler, geom, theta, amb);
theta = scatter_sysvec(theta, (K+H)\F);
theta=theta-273.15;% convert to degrees Celsius
% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',hot));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, T)));
geomT=nodal_field(struct ('name', ['geomT'], ...
    'data',[geom.values, theta.values]));
draw(fes1T3, gv, struct ('x',geomT, 'u',0*geomT,...
    'colorfield',colorfield, 'shrink',0.99));
draw(fes1T3, gv, struct ('x',geom, 'u',0*geom, ...
    'facecolor','none'));
draw(fes2T3, gv, struct ('x',geomT, 'u',0*geomT,...
    'colorfield',colorfield, 'shrink',0.99));
draw(fes2T3, gv, struct ('x',geom, 'u',0*geom, ...
    'facecolor','none'));
for iv=[-20:5:10]
draw_isosurface(fes1T3,gv, struct ('x', geomT,'u',0*geomT,'scalarfield',theta,'isovalue',iv,'color','k'));
draw_isosurface(fes2T3,gv, struct ('x', geomT,'u',0*geomT,'scalarfield',theta,'isovalue',iv,'color','k'));
end 
axis tight
set(gca,'DataAspectRatio', [1, 1, 100])
set_graphics_defaults(gcf)