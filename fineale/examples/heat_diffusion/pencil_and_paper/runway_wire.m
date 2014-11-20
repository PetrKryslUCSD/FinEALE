function runway_wire
% Runway with a heating wire.  Refined model.
u=physical_units_struct;
kappa=1.81*u.W/u.K/u.M; 
rho = 2350*u.KG/u.M^3;
cv =0.22*u.J/u.K/u.KG*rho;
d=  0.13*u.M;
% Tsoil=(8+  273.15)*u.K;
% Twire=(34+  273.15)*u.K;
% Ta=(-5+  273.15)*u.K;
Tsoil=(8)*u.K;
Twire=(34)*u.K;
Ta=(-5)*u.K;
h=15*u.W/u.K/u.M^2;
tolerance=d/10000;
num_integ_pts=1; % 1-point quadrature
[fens,fes, groups,edge_fes,edge_groups] = targe2_mesher({...
    ['curve 1 line '  num2str([0,0,d,0])],...
    ['curve 2 line '  num2str([d,0,d,d])],...
    ['curve 3 line '  num2str([d,d,d,2*d])],...
    ['curve 4 line '  num2str([d,2*d,0,2*d])],...
    ['curve 5 line '  num2str([0,2*d,0,d])],...
    ['curve 6 line '  num2str([0,d,0,0])],...
    ['subregion 1  property 1 ' ...
    '   boundary 1 2 3 4 5 6'],...
    ['m-ctl-point constant ' num2str(d/12)]
    }, 1.0);
prop = property_heat_diffusion(struct('thermal_conductivity',kappa,'source',0));
mater=material_heat_diffusion(struct('property',prop));
femm = femm_heat_diffusion(struct ('material',mater,...
    'fes',fes,...
    'integration_rule',tri_rule(struct('npts',num_integ_pts))));
sfemm = femm_heat_diffusion(struct ('material',mater,...
    'fes',subset(edge_fes,edge_groups{4}),...
    'integration_rule',gauss_rule(struct('dim',1,'order' ,2)),...
    'surface_transfer', h));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta=nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
fenids=connected_nodes  (subset(edge_fes,edge_groups{1}));
prescribed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+Tsoil;
theta = set_ebc(theta, fenids, prescribed, comp, val);
fenids=[intersect(connected_nodes(subset(edge_fes,edge_groups{2})),...
connected_nodes(subset(edge_fes,edge_groups{3})))];
prescribed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+Twire;%
theta = set_ebc(theta, fenids, prescribed, comp, val);
theta = apply_ebc (theta);
theta = numberdofs (theta);
amb = clone(theta, ['amb']);
fenids=connected_nodes  (subset(edge_fes,edge_groups{4}));
prescribed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+Ta;
amb = set_ebc(amb, fenids, prescribed, comp, val);
amb = apply_ebc (amb);


K = conductivity(femm, sysmat_assembler_sparse, geom, theta);
K = K + surface_transfer(sfemm, sysmat_assembler_sparse, geom, theta);
F = nz_ebc_loads_conductivity(femm, sysvec_assembler, geom, theta);
F = F + surface_transfer_loads(sfemm, sysvec_assembler, geom, theta, amb);
theta = scatter_sysvec(theta, K\F);

% Plot 
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
XY=geom.values;
conn= fes.conn;
draw(fes, gv, struct ('x',geom, 'u',0*geom, 'facecolor','none','edgecolor',[0.49, 0.49, 0.49]));
dcm=data_colormap(struct ('range',[min(T),max(T)], 'colormap',jet));
for  isovalue=linspace(min(T),max(T),20)
    draw_isosurface(femm,gv, struct ('x', geom,'u',0*geom,'scalarfield',theta,'isovalue',isovalue,'color',map_data(dcm, isovalue)));
end
draw_integration_points(femm, gv, struct ('x',geom,'u',0*geom, ...
    'theta', theta, 'scale', 00.000004,'color','red'));
draw_colorbar(gv, struct('colormap',dcm.colormap,...
                'position', [0.8, 0.1, 0.05, 0.6],...
                'minmax',[dcm.rmin,dcm.rmax],'label','Temperature'));
view(2)
set_graphics_defaults
labels('x','y')

end

