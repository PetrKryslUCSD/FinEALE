% Aluminum/steel assembly with rounded corners
%  ('Aluminum/steel assembly');
% Parameters:
Ea=70e9;
nua=0.33;
alphaa= 23e-6;
Es=200e9;
nus=0.3;
alphas= 12e-6;
thickness = 5.0;
integration_order =3;
scale = 0;
h=2.5;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    'curve 1 line 0 0 75 0',...
    'curve 2 line 75 0 75 40',...
    'curve 3 line 75 40 50 40',...
    'curve 4 line 50 40 50 17.5',...
    'curve 45 arc 47.5 15 50 17.5 center 47.5 17.5 rev',...
    'curve 5 line 47.5 15 0 15',...
    'curve 6 line 0 15 0 0',...
    'curve 7 line 50 40 0 40',...
    'curve 8 line 0 40 0 15',...
    ['subregion 1  property 1 boundary '...
    ' 1  2 3 4 -45 5 6'],...
    ['subregion 2  property 2 boundary '...
    ' -5 -4 7 8 45'],...
    ['m-ctl-point constant ' num2str(h)],...
    ['m-ctl-point 1 xy 50 15 near ' num2str(h/10)...
    ' influence  ' num2str(h/4)]
    }, thickness, struct('quadratic',true));
% Material
propa = property_deformation_linear_iso ...
    (struct('E',Ea,'nu',nua,'alpha', alphaa));
matera = material_deformation_linear_biax (struct('property',propa, ...
    'reduction','stress'));
props = property_deformation_linear_iso ...
    (struct('E',Es,'nu',nus,'alpha', alphas));
maters = material_deformation_linear_biax (struct('property',props, ...
    'reduction','stress'));
% Finite element block
femma = femm_deformation_linear (struct ('material',matera,...
    'fes',subset(fes,groups{2}),...
    'integration_rule',tri_rule (struct('npts',3))));
femms = femm_deformation_linear (struct ('material',maters,...
    'fes',subset(fes,groups{1}),...
    'integration_rule',tri_rule (struct('npts',3))));
geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[0 0 0 40],'inflate', 0.01));
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=1*ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
ebc_fenids=fenode_select (fens,struct('box',[0 75 40 40],'inflate', 0.01));
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=2*ones(1,length (ebc_fenids));
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
u   = apply_ebc (u);
u   = numberdofs (u);
% Temperature field
dT = nodal_field(struct ('name',['dT'], 'dim', 1, ...
    'data',zeros(count(fens),1)+70));
% Assemble the system matrix
K = stiffness(femms, sysmat_assembler_sparse,  geom, u)...
    +stiffness(femma, sysmat_assembler_sparse, geom, u);
% Load
F = thermal_strain_loads(femma, sysvec_assembler, geom, u, dT)...
    +thermal_strain_loads(femms, sysvec_assembler, geom, u, dT);
% Solve
u = scatter_sysvec(u, K\F);
% get(u,'values')

% Plot
gv=graphic_viewer;
gv=reset (gv,struct ('limits', [0, 100, -8, 40]));
set(gca,'FontSize', 12)
cmap=jet;
cmpn=3;
flda = field_from_integration_points(femma, geom, u, dT, 'Cauchy', cmpn);
flds = field_from_integration_points(femms, geom, u, dT, 'Cauchy', cmpn);
nvalsa=flda.values;
nvalss=flds.values;
nvalmin =min(min(nvalsa),min(nvalss));
nvalmax =max(max(nvalsa),max(nvalss))
dcm=data_colormap(struct ('range',[nvalmin,nvalmax], 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvalsa)));
draw(femma, gv, struct ('x',geom,'u', scale*u, 'colorfield',colorfield, 'shrink',1));
% draw(femma, gv, struct ('x',geom,'u', 0*u, 'facecolor',' none'));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvalss)));
draw(femms, gv, struct ('x',geom,'u', scale*u, 'colorfield',colorfield, 'shrink',1));
% draw(femms, gv, struct ('x',geom,'u', 0*u, 'facecolor',' none'));
lighting  none;
title (['h=' num2str(h)]);
colormap(cmap);
draw_colorbar(gv,struct('position',[0.72 0.33 0.05 0.5],...
    'minmax',[nvalmin,nvalmax],'label','\sigma_{xy}'));
view (2)
%  saveas(gcf, [mfilename '-' num2str(h) '.png'], 'png');
    

