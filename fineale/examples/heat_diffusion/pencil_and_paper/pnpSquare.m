% Square domain with internal heat source  and fixed temperature 
% on the boundary.
function pnpSquare
kappa=[1.8 0; 0 1.8]; % conductivity matrix
Q=15; % uniform heat source
a=1;
num_integ_pts=1; % 1-point quadrature
[fens,fes] = targe2_mesher({... 
        ['curve 1 line ' num2str(-2*a) ' 0 0 0'],...
        ['curve 2 line 0 0 ' num2str(-2*a) ' ' num2str(2*a) ''],...
        ['curve 3 line ' num2str(-2*a) ' ' num2str(2*a) ' ' num2str(-2*a) ' 0'],...
        'subregion 1  property 1 boundary 1 2 3 ',...
        ['m-ctl-point constant 0.12']
        }, 1.0);
prop=property_heat_diffusion(struct('thermal_conductivity',kappa,'source',Q));
mater=material_heat_diffusion (struct('property',prop));
femm = femm_heat_diffusion (struct ('material',mater,...
    'fes',fes,...
    'integration_rule',tri_rule(struct('npts', num_integ_pts))));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta = nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
fenids=[fenode_select(fens,struct('box',[-2 0 0 0],...
    'inflate', 0.01))];
fixed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+20;
theta = set_ebc(theta, fenids, fixed, comp, val);
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
axis equal vis3d
set_graphics_defaults(gcf)
(theta.values-20)/(Q*a^2/6/kappa(1,1));