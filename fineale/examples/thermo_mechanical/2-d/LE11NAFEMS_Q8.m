% NAFEMS LE11 benchmark with Q8 elements.
% This is a test recommended by the National Agency for Finite Element
% Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
% Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.  
%
% Target solution: Direct stress,  = –105 MPa at point A.
function  LE11NAFEMS_Q8
% Parameters:
pu= physical_units_struct;
Ea= 210000*pu.MEGA*pu.PA;
nua= 0.3;
alphaa=2.3e-4;
sigmaA=-105*pu.MEGA*pu.PA;
nref= 0;
    X=[1.    , 0.;%A
 1.4   , 0.;%B
 0.995184726672197   0.098017140329561;
 1.393258617341076 0.137223996461385;
 0.980785,0.195090;%
 1.37309939,0.27312645;
  0.956940335732209   0.290284677254462
  1.339716470025092 0.406398548156247
 0.9238795, 0.38268;%C
 1.2124, 0.7;%D
 0.7071, 0.7071;%E
 1.1062, 1.045;%F
  0.7071, (0.7071+1.79)/2;%(E+H)/2
1.    , 1.39;%G
 0.7071, 1.79;%H
 1.    , 1.79;%I
]*pu.M;
tolerance =1e-6;
fens=fenode_set(struct('xyz',X));
fes=fe_set_Q4(struct('conn',[1,2,4,3;3,4,6,5;5,6,8,7;7,8,10,9;9,10,12,11;11,12,14,13;13,14,16,15], 'axisymm', true));
for ref=1:nref
    [fens,fes]=Q4_refine(fens,fes);
    list=fenode_select(fens,struct('distance',1.0+0.1/2^nref, 'from',[0,0],'inflate', tolerance));
    fens= onto_sphere(fens,1.0,list);
end
[fens,fes]=Q4_to_Q8(fens,fes,struct('axisymm', true));
list=fenode_select(fens,struct('distance',1.0+0.1/2^nref, 'from',[0,0],'inflate', tolerance));
fens= onto_sphere(fens,1.0,list);
 %     drawmesh({fens,fes},'fes','facecolor','red');
 %     view(2)
% Material
propa = property_deformation_linear_iso ...
    (struct('E',Ea,'nu',nua,'alpha', alphaa));
matera = material_deformation_linear_biax (struct('property',propa, ...
    'reduction','axisymm'));

% Finite element block
femma = femm_deformation_linear (struct ('material',matera,...
    'fes',fes,...
    'integration_rule',gauss_rule (struct('dim',2, 'order',3))));

geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=[fenode_select(fens,struct('box',[0 1.4 0 0],'inflate', tolerance)),...
    fenode_select(fens,struct('box',[0 1.4 1.79,   1.79],'inflate', tolerance))];
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=2*ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);

u   = apply_ebc (u);
u   = numberdofs (u);
% Temperature field
x=fens.xyz;
dT = nodal_field(struct ('name',['dT'], 'dim', 1, ...
    'data',x(:,1)+x(:,2)));
% Assemble the system matrix
K = stiffness(femma, sysmat_assembler_sparse, geom, u);
% Load
F = thermal_strain_loads(femma, sysvec_assembler, geom, u, dT);
% Solve
u = scatter_sysvec(u, K\F);
% get(u,'values')

% Plot
gv=graphic_viewer;
gv=reset (gv,struct ('limits', [0, 2, 0, 1.8]));
set(gca,'FontSize', 12)
cmap=jet;
cmpn=2;
flda = field_from_integration_points(femma, geom, u, dT, 'Cauchy', cmpn);
nvalsa=flda.values/(pu.MEGA*pu.PA);
nvalmin =min(nvalsa);
nvalmax =max(nvalsa);
dcm=data_colormap(struct ('range',[nvalmin,nvalmax], 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvalsa)));
draw(femma, gv, struct ('x',geom,'u', 100*u, 'colorfield',colorfield, 'shrink',1));
draw(mesh_boundary(femma.fes,struct('axisymm', true,'other_dimension', 0.1)), gv, struct ('x',geom,'u', 0*u, 'edgecolor','r'));
% draw(femms, gv, struct ('x',geom,'u', 0*u, 'facecolor',' none'));
lighting  none;
colormap(cmap);
draw_colorbar(gv,struct('position',[0.72 0.33 0.05 0.5],...
    'minmax',[nvalmin,nvalmax],'label','\sigma_{y}'));
view (2)
%  saveas(gcf, [mfilename '-' num2str(h) '.png'], 'png');
    
nA =fenode_select(fens,struct('box',[0 1 0 0],'inflate', tolerance));
disp(['Stress at point A: ' num2str(nvalsa(nA)) ', i. e.  ' num2str(nvalsa(nA)*pu.MEGA*pu.PA/sigmaA*100) ' %'])
end