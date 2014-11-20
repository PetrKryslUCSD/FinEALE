% Square domain with a square inclusion of orthotropic material.
function pnpSquareInSquare
N = 12;
N_f = 4;
node2dof=zeros(1,N);
node2dof([9,7,4,6]) =1:N_f;
conninner = [9,7,6; 7,4,6];
connouter = [1,4,2; 4,3,2; 4,7,3; 1,6,4; 1,5,6; 5,10,6; 10,9,6;...
    10,11,9; 11,12,9; 12,7,9; 12,8,7; 8,3,7];
x= [-48,48; 0,48; 48,48; 0,31; -48,0;...
    -31,0; 31,0; 48,0; 0,-31; -48,-48; ...
    0,-48; 48,-48];
T20=20; T57 =57;% Boundary conditions, degrees Celsius
pT =sym(zeros(N,1));
pT([10, 11, 12]) =T20;% bottom horizontal face
pT([1, 2, 3]) =T57;% top horizontal face

kappainner=[2.25 0; 0 0.06]; % orthotropic conductivity matrix
kappaouter=[0.25 0; 0 0.25]; % isotropic conductivity matrix
alpha =-45;% local material orientation angle
ca=cos(2*pi/360*alpha); sa=sin(2*pi/360*alpha);
Rm = [ca, -sa;sa, ca];% local material directions

fens =fenode_set(struct('xyz',x));
fesinner=fe_set_T3(struct('conn',conninner));
fesouter=fe_set_T3(struct('conn',connouter));

propinner = property_heat_diffusion(struct('thermal_conductivity',kappainner,...
    'source',0));
materinner=material_heat_diffusion(struct('property',propinner));
femminner = femm_heat_diffusion(struct ('material',materinner,...
    'fes',fesinner,...
    'integration_rule',tri_rule(struct('npts', 1)),'Rm',Rm));
propouter = property_heat_diffusion(struct('thermal_conductivity',kappaouter,...
    'source',0));
materouter=material_heat_diffusion(struct('property',propouter));
femmouter = femm_heat_diffusion(struct ('material',materouter,...
    'fes',fesouter,...
    'integration_rule',tri_rule(struct('npts', 1))));
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
theta=nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
fenids=[fenode_select(fens,struct('box',[-48 48 -48 -48],...
    'inflate', 0.01))];
fixed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+20;% ambient temperature
theta = set_ebc(theta, fenids, fixed, comp, val);
fenids=[fenode_select(fens,struct('box',[-48 48 48 48],...
    'inflate', 0.01))];
fixed=ones(length(fenids),1);
comp=[];
val=zeros(length(fenids),1)+57;% hot inner surface
theta = set_ebc(theta, fenids, fixed, comp, val);
theta = apply_ebc (theta);
theta = numberdofs (theta);
K = conductivity(femminner, sysmat_assembler_sparse, geom, theta)+ ...
    conductivity(femmouter, sysmat_assembler_sparse, geom, theta);
F = nz_ebc_loads_conductivity(femminner, sysvec_assembler, geom, theta)+ ...
    nz_ebc_loads_conductivity(femmouter, sysvec_assembler, geom, theta);
theta = scatter_sysvec(theta, K\F);

% Plothot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
dcm=data_colormap(struct ('range',[min(T),max(T)], 'colormap',hot(9)));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, T)));
geomT=nodal_field(struct ('name', ['geomT'], ...
    'data',[geom.values,theta.values]));
draw(femmouter, gv, struct ('x',geomT, 'u',0*geomT,...
    'colorfield',colorfield, 'shrink',1));
draw(femmouter, gv, struct ('x',geom, 'u',0*geom, 'facecolor','none'));
draw(femminner, gv, struct ('x',geomT, 'u',0*geomT,...
    'colorfield',colorfield, 'shrink',1));
draw(femminner, gv, struct ('x',geom, 'u',0*geom, 'facecolor','none'));
for  isovalue=linspace(min(T),max(T),20)
    draw_isosurface(femminner,gv, struct ('x', geomT,'u',0*geomT,'scalarfield',theta,'isovalue',isovalue,'color','k'));
    draw_isosurface(femmouter,gv, struct ('x', geomT,'u',0*geomT,'scalarfield',theta,'isovalue',isovalue,'color','k'));
end
axis equal
xlabel x, ylabel y, zlabel T 
set_graphics_defaults
view(2)