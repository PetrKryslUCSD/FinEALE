% Thermomechanical analysis of thermal-expansion driven MEMS actuator
actuator2;% run the thermal analysis
E= 169000000;
nu =0.278;
alpha = 2.5e-6;
% Material: mechanical
mechprop = property_deformation_linear_iso (struct ('E',E,...
    'nu',nu,'alpha', alpha));
mechmater=material_deformation_linear_triax(struct('property',mechprop));

% Finite element block: mechanical
mechfemm = femm_deformation_linear (struct ('material',mechmater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct( 'dim',3,'order',2))));

u=0*geom;
fenids=fenode_select(fens,struct('box',[x0,x4,y0,y0,z0,z3],...
    'inflate', t/1000)) ; % fixed displacement
fixed=fenids*0+1;
comp=[];
val=zeros(length(fenids),1)+0;
u = set_ebc(u, fenids, fixed, comp, val);
fenids=fenode_select(fens,struct('box',[x0,x0,y0,100*y3,z0,100*z3],...
    'inflate', t/1000)) ; % fixed displacement
fixed=fenids*0+1;
comp=fenids*0+1;% x-compe
val=zeros(length(fenids),1)+0;
u = set_ebc(u, fenids, fixed, comp, val);
u = apply_ebc (u);

u = numberdofs (u);
K = stiffness(mechfemm, sysmat_assembler_sparse, geom, u);
F =  thermal_strain_loads(mechfemm, sysvec_assembler, geom, u, theta-T_substrate);
u = scatter_sysvec(u, K\F);

% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
draw(mechfemm, gv, struct ('x',geom, 'u',0*u,...
    'facecolor','none'));
U=u.values;
Uz=U(:,3);
cmap= jet;
dcm=data_colormap(struct('range',[min(Uz),max(Uz)],'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, Uz)));
draw(mechfemm, gv, struct ('x',geom, 'u',u,...
    'colorfield',colorfield));
draw_colorbar(gv, struct('colormap',cmap,'label','Uz',...
    'minmax',[min(Uz),max(Uz)]))

interact(gv);
