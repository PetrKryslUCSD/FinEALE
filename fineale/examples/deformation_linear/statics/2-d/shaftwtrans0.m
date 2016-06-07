% Axially symmetric analysis of a shaft: axis(1) + clamped face(2) 
disp('Axially symmetric analysis of a shaft: axis(1) + clamped face(2)');


% Parameters:
E=210e6;% MPa
nu=0.29;
% geometry
R1= 20;% millimeters
Rf=5;% millimeters
R2=R1-Rf;
L1=40;%mm
L2=L1;%mm
% fixed force
forc=pi*R1^2*200;%N

% Mesh'
mesh_size=3;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    ['curve 1 line 0 0 ' num2str(R1) ' 0'],...
    ['curve 2 line ' num2str(R1) ' 0 ' num2str(R1) ' ' num2str(L1) ],...
    ['curve 3 arc ' num2str(R1) ' ' num2str(L1) ' ' num2str(R2) ' ' num2str(L1+Rf) ' center ' num2str(R1) ' ' num2str(L1+Rf)],...
    ['curve 4 line ' num2str(R2) ' ' num2str(L1+Rf) '  ' num2str(R2) ' ' num2str(L1+L2)],...
    ['curve 5 line ' num2str(R2) ' ' num2str(L1+L2) ' ' num2str(0) ' ' num2str(L1+L2)],...
    ['curve 6 line '  num2str(0) ' ' num2str(L1+L2)  ' ' num2str(0) ' ' num2str(0)],...
    ['subregion 1  property 1 boundary 1 2 3 4 5 6'],...
    ['m-ctl-point constant ' num2str(mesh_size)],...
    ['m-ctl-point 1 xy ' num2str(R1) ' ' num2str(L1+Rf) ' near ' num2str(mesh_size/200) ' influence ' num2str(mesh_size/20)],...
    ['m-ctl-point 1 xy ' num2str(R1) ' ' num2str(0) ' near ' num2str(mesh_size/300) ' influence ' num2str(mesh_size/60)],...
    }, 1.0, struct('axisymm',true,'quadratic',true));
% drawmesh({fens,fes},'fes','facecolor','red')

% Material
prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_linear_biax (struct('property',prop, ...
    'reduction','axisymm'));
% Finite element block
femm = femm_deformation_linear(struct ('material',mater, 'fes',fes,...
    'integration_rule',tri_rule (struct('npts',3))));
efemm = femm_deformation_linear (struct ('material',mater, 'fes',subset(edge_fes,edge_groups{5}),...
    'integration_rule',gauss_rule (struct('dim',1,'order', 2))));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[0 R1 0 0],'inflate',R1/10000));
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=[];
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
% The axis of symmetry
ebc_fenids=fenode_select (fens,struct('box',[0 0 0 L1+L2],'inflate',L1/10000));
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=ebc_fixed*0+1;
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);
% Assemble the system matrix
K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
% Load
F = nz_ebc_loads(femm, sysvec_assembler, geom, u);
fi=force_intensity(struct('magn',[0;forc/(pi*R2^2)]));
F = F + distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
% Solve
u = scatter_sysvec(u, K\F);

% Plot
gv=graphic_viewer;
gv=reset (gv,struct ('limits',[0 1.06*R1 0 1.1*(L1+L2)]));
scale=1;
cmap = jet;
fld = field_from_integration_points(femm, geom, u, [], [], [], 'Cauchy',2);
nvals=fld.values;%min(nvals),max(nvals)
nvalsrange=[min(nvals),max(nvals)];
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
geom3=nodal_field(struct ('name', ['geom3'], ...
    'data',[geom.values,0.25*nvals]));
u3=nodal_field(struct ('name', ['u3'], ...
    'data',[u.values,0*nvals]));
draw(femm,gv, struct ('x', geom3, 'u', +scale*u3,'colorfield',colorfield, 'shrink',1.0));
draw_colorbar(gv, struct('colormap',cmap,'position',[0.8 0.15 0.05 0.7],...
    'minmax',nvalsrange,'label','\sigma_y'));
% view (2)
