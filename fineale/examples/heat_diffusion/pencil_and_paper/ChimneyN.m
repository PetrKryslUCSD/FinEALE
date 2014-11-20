u= physical_units_machine;
kappa=0.5*eye(2)*u('W/m/K'); % conductivity matrix
hi=5*u('W/m^2/K');
ho=15*u('W/m^2/K');
a=0.2*u('m');
Ti = 1000*u('K');
To = -20*u('K');
fens=fenode_set(struct('id',(1:4)','xyz',[0,0; 2*a,0; 0,a; a,a]));
fes=fe_set_T3(struct('id',(1:2)','conn',[3,1,4; 4,1,2]));
edge_fes=fe_set_L2(struct('id',(1:2)','conn',[1,2; 4,3],'other_dimension', 1));

prop = property_heat_diffusion(struct('thermal_conductivity',kappa,'source',0));
mater=material_heat_diffusion(struct('property',prop));
femm = femm_heat_diffusion(struct ('material',mater,...
    'fes',fes,...
    'integration_rule',tri_rule(struct('npts',1))));
iedgefemm = femm_heat_diffusion (struct ('mater',mater,...
    'fes',subset(edge_fes,2),...
    'integration_rule',gauss_rule(struct('dim',1,'order' ,2)),...
    'surface_transfer', hi));
oedgefemm = femm_heat_diffusion (struct ('mater',mater,...
    'fes',subset(edge_fes,1),...
    'integration_rule',gauss_rule(struct('dim',1,'order' ,2)),...
    'surface_transfer', ho));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta=nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',...
    geom.nfens));
fenids=1:count(fens);
theta = numberdofs (theta);
iamb = clone(theta, ['iamb']);
prescribed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+Ti;
iamb = set_ebc(iamb, fenids, prescribed, comp, val);
iamb = apply_ebc (iamb);
oamb = clone(theta, ['oamb']);
prescribed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+To;
oamb = set_ebc(oamb, fenids, prescribed, comp, val);
oamb = apply_ebc (oamb);
K =  conductivity(femm, sysmat_assembler_sparse, geom, theta);
full (K)
H = surface_transfer(iedgefemm, sysmat_assembler_sparse, geom, theta)+...
     surface_transfer(oedgefemm, sysmat_assembler_sparse, geom, theta);
full(H)
F = surface_transfer_loads(iedgefemm, sysvec_assembler, geom, theta, iamb)+...
     surface_transfer_loads(oedgefemm, sysvec_assembler, geom, theta, oamb);
 F
theta = scatter_sysvec(theta, (K+H)\F);
T = theta.values


% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',hot));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, T)));
geomT=nodal_field(struct ('name', ['geomT'], ...
    'data',[geom.values, theta.values/100]));
draw(fes, gv, struct ('x',geomT, 'u',0*geomT,...
    'colorfield',colorfield, 'shrink',01));
draw(fes, gv, struct ('x',geom, 'u',0*geom, ...
        'facecolor','none'));

axis equal vis3d
zlabel('T/100')
format long

% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',hot));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, T)));
for  isovalue=[0, 100, 200, 300, 400, 500, 600]
    draw_isosurface(femm,gv, struct ('x', geom,'u',0*geom,'scalarfield',theta,'isovalue',isovalue,'color',map_data(dcm, isovalue)));
end
draw(fes, gv, struct ('x',geom, 'u',0*geom, ...
        'facecolor','none'));
draw_integration_points(femm, gv, struct ('x',geom,'u',0*geom, ...
    'theta', theta, 'scale', 00.0002,'color','red'));
draw_colorbar(gv, struct('colormap',dcm.colormap,...
                'position', [0.8, 0.1, 0.05, 0.6],...
                'minmax',[dcm.rmin,dcm.rmax],'label','Temperature'));
view (2)
axis equal 
zlabel('T/100')
format long


