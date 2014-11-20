%% Solid cylinder/taper/sphere—-temperature loading
%

%%
% Link to the  <matlab:edit('pub_LE11NAFEMS_Q8') m-file>.
%

%% Description
%
% The solid cylinder/taper/sphere axially-symmetric part represented in
% Figure 1 is exposed to linearly varying temperature in the plane of the
% cross-section. The temperature in the coordinates $r$  (the radial
% coordinate) and $z$ (the axial ccoordinate)  is given as $T=r+z$. The
% goal is to find  the mechanical stress at the point A induced by the
% thermal expansion.
%

%%
%
% <html> <table border=0><tr><td> <img
% src="../docs/pub_LE11NAFEMS_Q8_2.jpg"> </td></tr> <tr><td>Figure 1.
% Definition of the geometry of the part</td></tr> </table> </html>

%%
% The part is constrained against axial expansion along the lines of HI and
% AB. The Young's modulus is 210 GPa, the Poisson's ratio is .3, and the
% coefficient of thermal expansion is 2.3e-4/degree Celsius.
%%
%
% This is a test recommended by the National Agency for Finite Element
% Methods and Standards (U.K.): Test LE11 from NAFEMS Publication TNSB,
% Rev. 3, “The Standard NAFEMS Benchmarks,” October 1990.
%
% Target solution: Compressive  axial stress $\sigma_z$  = –105 MPa at
% point A.



%% Solution
%
function  pub_LE11NAFEMS_Q8
    %%
    % The toolkit has a helpful physical-units facility.  The variable |pu|
    % is a structure with fields that define basic  units and basic
    % multipliers (for instance, mega).
    pu= physical_units_struct;
    %%
    % Set the material properties.
    Ea= 210000*pu.MEGA*pu.PA;% Young's modulus
    nua= 0.3;% Poisson ratio
    alphaa=2.3e-4;% coefficient of thermal expansion
    
    %%
    % This is the target stress value.
    sigma_z_A_ref=-105*pu.MEGA*pu.PA;
    
    %%
    % The mesh  will be created in a very coarse representation from the
    % key points in the drawing.
    rz=[1.    , 0.;%A
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
    fens=fenode_set(struct('xyz',rz));
    fes=fe_set_Q4(struct('conn',[1,2,4,3;3,4,6,5;5,6,8,7;7,8,10,9;9,10,12,11;11,12,14,13;13,14,16,15], 'axisymm', true));
    
    %%
    % If needed, the initial mesh  can be refined by bisection.  Just set
    % |nref| greater than zero.  Note that  the nodes located along the
    % edges are moved onto the  spherical surface when they _should be_ on
    % the spherical surface.  This is important in order to ensure
    % convergence to the proper value of the stress.  Just refining  the
    % initial mesh without repositioning of the nodes would mean that the
    % refinement would preserve a concave corner where in reality there is
    % none.  The stress would be artificially raised and convergence would
    % not be guaranteed.
    
    nref= 0;
    for ref=1:nref
        [fens,fes]=Q4_refine(fens,fes);
        list=fenode_select(fens,struct('distance',1.0+0.1/2^nref, 'from',[0,0],'inflate', tolerance));
        fens= onto_sphere(fens,1.0,list);
    end
    
    %%
    % The mesh is now converted to the serendipity eight-node elements.
    % Note that their attribute |axisymm| is set  to reflect  the axial
    % symmetry of the problem.
    [fens,fes]=Q4_to_Q8(fens,fes,struct('axisymm', true));
    
    %%
    % The nodes within the radial distance of 1.0 of the origin (i. e.
    % tthose on the spherical surface)  are repositioned one more time to
    % be located on the spherical surface for sure. (Recall  that we have
    % inserted additional nodes at the midpoints of the edges when the mesh
    % was converted to quadratic elements.)
    list=fenode_select(fens,struct('distance',1.0+0.1/2^nref, 'from',[0,0],'inflate', tolerance));
    fens= onto_sphere(fens,1.0,list);
    
    %%
    % The mesh is drawn as a check.
    drawmesh({fens,fes},'fes','facecolor','red');
    view(2)
    
    
    %%
    % We are ready to create the  finite element model machine and to use
    % it to construct  the global system for the displacements.
    %%
    % The material is created from the property object.  Note that the
    % |alpha| attribute is the thermal expansion coefficient.
    propa = property_deformation_linear_iso ...
        (struct('E',Ea,'nu',nua,'alpha', alphaa));
    
    %%
    % The  model-order reduction that accounts for the axial symmetry needs
    % to be reflected in the mmaterial stiffness matrix  that's the
    % material object computes.  The |reduction| property needs to be set.
    matera = material_deformation_linear_biax (struct('property',propa, ...
        'reduction','axisymm'));
    
    
    %%
    % The finite element  model machine puts together the material, the
    % finite elements,  and the integration rule. The Gauss quadrature
    % with 3 x 3 points  gives good accuracy in this case. Compared with 2
    % x 2 quadrature to appreciate the difference.
    femma = femm_deformation_linear (struct ('material',matera,...
        'fes',fes,...
        'integration_rule',gauss_rule (struct('dim',2, 'order',3))));
    
    %%
    % The geometry nodal field is created from the node set.   The
    % displacement field is created by cloning the geometry and then
    % zeroing out the nodal parameters.
    geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
    u   = clone(geom,'u');
    u   = u*0; % zero out
    
    %%
    % The EBCs are applied the next.  Only the axial degrees of freedom at
    % the bottom and top are fixed to zero.
    ebc_fenids=[fenode_select(fens,struct('box',[0 1.4 0 0],'inflate', tolerance)),...
        fenode_select(fens,struct('box',[0 1.4 1.79,   1.79],'inflate', tolerance))];
    ebc_fixed=ones(1,length (ebc_fenids));% The degrees of freedom are being fixed.
    ebc_comp=2*ones(1,length (ebc_fenids));% The axial component is being fixed
    ebc_val=ebc_fenids*0;% Set to zero for the prescribed degrees of freedom.
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    
    %%
    % The EBCs are applied and the degrees of freedom are numbered.
    u   = apply_ebc (u);
    u   = numberdofs (u);
    
    %%
    % We create the temperature field using the formula.
    dT = nodal_field(struct ('name',['dT'], 'dim', 1, ...
        'data',fens.xyz(:,1)+fens.xyz(:,2)));
    
    %%
    % And we are ready to assemble the system matrix.
    K = stiffness(femma, sysmat_assembler_sparse, geom, u);
    
    %%
    % The mechanical loads are computed from the thermal strains.
    F = thermal_strain_loads(femma, sysvec_assembler, geom, u, dT);
    
    %%
    % And  the solution for the free degrees of freedom is obtained.
    u = scatter_sysvec(u, K\F);
    
    
    %%
    % The stress  is recovered from the stress calculated at the
    % integration points.  The method |field_from_integration_points| uses
    % inverse-distance interpolation to compute the nodal stress from the
    % quadrature-point stresses.
    cmpn=2;% this is the axial stress
    flda = field_from_integration_points(femma, geom, u, dT, 'Cauchy', cmpn);
    
    %%
    % Now that we have the nodal field  for the axial stress, we can plot
    % the axial stress painted on the deformed geometry.
    figure('visible', 'off')
    gv=graphic_viewer;
    gv=reset (gv,struct ('limits', [0, 2, 0, 1.8]));
    nvalsa=flda.values/(pu.MEGA*pu.PA);
    nvalmin =min(nvalsa);
    nvalmax =max(nvalsa);
    dcm=data_colormap(struct ('range',[nvalmin,nvalmax], 'colormap',jet));
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvalsa)));
    draw(femma, gv, struct ('x',geom,'u', 100*u, 'colorfield',colorfield, 'shrink',1));
    %%
    % The undeformed geometry is shown using an outline.
    draw(mesh_boundary(femma.fes,struct('axisymm', true,'other_dimension', 0.1)), gv, struct ('x',geom,'u', 0*u, 'edgecolor','r'));
    
    draw_colorbar(gv,struct('position',[0.72 0.33 0.05 0.5],...
        'minmax',[nvalmin,nvalmax],'label','\sigma_{y}'));
    view (2)
    set(gcf,'visible', 'on')
    %%
    % The  computed stress at the node that is located at the point A  is
    % going to be now extracted from the nodal field for the stress.
    nA =fenode_select(fens,struct('box',[0 1 0 0],'inflate', tolerance));
    disp(['Stress at point A: ' num2str(nvalsa(nA)) ', i. e.  ' num2str(nvalsa(nA)*pu.MEGA*pu.PA/sigma_z_A_ref*100) ' %'])
    
    %% Discussion
    %
    %%
    % The accuracy of the computed stress depends on the fineness  of the
    % mesh and also on the number of quadrature points.
    %%
    % # The fineness of the mesh affects the general accuracy of the
    % solution.
    % # The number of quadrature points affects the extrapolation procedure
    % from the quadrature points to the nodes to create the smooth stress
    % field.
end