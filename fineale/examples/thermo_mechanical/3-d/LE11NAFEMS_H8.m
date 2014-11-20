% NAFEMS LE11 benchmark, 3-D hexahedral mesh
% This is a test recommended by the National Agency for Finite Element
% Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
% Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
%
% Target solution: Direct stress,  = –105 MPa at point A.
function  LE11NAFEMS
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
    tolerance =1e-3;
    fens=fenode_set(struct('xyz',X));
    fes=fe_set_Q4(struct('conn',[1,2,4,3;3,4,6,5;5,6,8,7;7,8,10,9;9,10,12,11;11,12,14,13;13,14,16,15], 'axisymm', true));
    for ref=1:nref
        [fens,fes]=Q4_refine(fens,fes);
        list=fenode_select(fens,struct('distance',1.0+0.1/2^nref, 'from',[0,0],'inflate', tolerance));
        fens= onto_sphere(fens,1.0,list);
    end
    nLayers=1;
    angslice =pi/16;
    [fens,fes] = H8_extrude_Q4(fens,fes,nLayers,@(xyz,k)[xyz,0]-(k-1/2)/nLayers*[0,0,angslice]);
    bfes=mesh_boundary(fes);
    f1l=fe_select (fens,bfes,struct('box',[-inf,inf,-inf,inf,-angslice/2,-angslice/2],'inflate',100*eps));
    f2l=fe_select (fens,bfes,struct('box',[-inf,inf,-inf,inf,angslice/2,angslice/2],'inflate',100*eps));
    fens = transform_apply(fens,@(x,d)[x(1)*cos(x(3)),x(2),-x(1)*sin(x(3))], []);
    gv=drawmesh({fens,fes},'fes','facecolor','none');
    gv=drawmesh({fens,subset(bfes,f1l)},'gv',gv,'fes','facecolor','y');
    gv=drawmesh({fens,subset(bfes,f2l)},'gv',gv,'fes','facecolor','r');
    labels  ([])
    %     view(2)
    % Material
    propa = property_deformation_linear_iso ...
        (struct('E',Ea,'nu',nua,'alpha', alphaa));
    matera = material_deformation_linear_triax (struct('property',propa));
    
    % Finite element block
    femma = femm_deformation_linear (struct ('material',matera,...
        'fes',fes,...
        'integration_rule',gauss_rule (struct('dim',3, 'order',2))));
    sfemma = femm_deformation_linear (struct ('material',matera,...
        'fes',subset(bfes,[f1l,f2l]),...
        'integration_rule',gauss_rule (struct('dim',2, 'order',2)),...
        'surface_normal_spring_coefficient',1000*Ea));
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    u   = clone(geom,'u');
    u   = u*0; % zero out
    % Apply EBC's
    ebc_fenids=[fenode_select(fens,struct('box',[-inf,inf, 0,0,-inf,inf,],'inflate', tolerance)),...
        fenode_select(fens,struct('box',[-inf,inf, 1.79,1.79,-inf,inf,],'inflate', tolerance))];
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=2*ones(1,length (ebc_fenids));
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    
    u   = apply_ebc (u);
    u   = numberdofs (u);
    % Temperature field
    x=fens.xyz;
    dT = nodal_field(struct ('name',['dT'], 'dim', 1, ...
        'data',sqrt(x(:,1).^2+x(:,3).^2)+x(:,2)));
    % Assemble the system matrix
    K = stiffness(femma, sysmat_assembler_sparse, geom, u);
    H = surface_normal_spring_stiffness(sfemma, sysmat_assembler_sparse, geom, u);
    % Load
    F = thermal_strain_loads(femma, sysvec_assembler, geom, u, dT);
    % Solve
    u = scatter_sysvec(u, (K+H)\F);
    % get(u,'values')
    
    % Plot
    gv=graphic_viewer;
    gv=reset (gv,struct ('limits', [0, 2, 0, 1.8]));
    set(gca,'FontSize', 12)
    cmap=jet;
    cmpn=2;
    % flda = dT; nvalsa=flda.values;
    flda = field_from_integration_points(femma, geom, u, dT, 'Cauchy', cmpn);
    nvalsa=flda.values/(pu.MEGA*pu.PA);
    nvalmin =min(nvalsa);
    nvalmax =max(nvalsa);
    dcm=data_colormap(struct ('range',[nvalmin,nvalmax], 'colormap',cmap));
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvalsa)));
    draw(femma, gv, struct ('x',geom,'u', 100*u, 'colorfield',colorfield, 'shrink',1));
    % draw(mesh_boundary(femma.fes,struct('axisymm',
    % true,'other_dimension', 0.1)), gv, struct ('x',geom,'u', 0*u,
    % 'edgecolor','r'));
    draw(femma, gv, struct ('x',geom,'u', 0*u, 'facecolor',' none'));
    lighting  none;
    colormap(cmap);
    draw_colorbar(gv,struct('position',[0.72 0.33 0.05 0.5],...
        'minmax',[nvalmin,nvalmax],'label','\sigma_{y}'));
    labels  ([])
    %  saveas(gcf, [mfilename '-' num2str(h) '.png'], 'png');
    
    nA =intersect (fenode_select(fens,struct('box',[-inf,inf, 0 0 -inf,inf],'inflate', tolerance)),...
        fenode_select(fens,struct('distance',1.0+0.01, 'from',[0,0,0],'inflate', tolerance)));
    sigAh=mean(nvalsa(nA));
    disp(['Stress at point A: ' num2str(sigAh) ', i. e.  ' num2str(sigAh*pu.MEGA*pu.PA/sigmaA*100) ' %'])
end