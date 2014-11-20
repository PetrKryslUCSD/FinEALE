%% Heated brick.  Solution with  quadratic tetrahedra.
%

%%
% Link to the  <matlab:edit('pub_heated_brick_tet') m-file>.
%

%% Description
%
% The brick represented in Figure 1 has insulated lateral surfaces (BCGF,
% CDHG, DAEH, ABFE).  The  surface  ABCD  is  maintained  at  constant
% temperature  T  =  0  degrees,  and  the temperature is linearly
% distributed at EFGH surface so that T = 0 degrees at E point, T = 1
% degrees at H and T = 2 degrees at point F. Find temperature at the point
% P. The reference solution for the temperature at point P is 0.901 degrees. Thermal
% conductivity $\kappa = 1$.
%

%%
% The problem was solved with hexahedra in tutorials pub_heated_brick and
% pub_heated_brick_quadratic. In this case quadratic tetrahedra will be
% used.
%%
%
% <html>
% <table border=0><tr><td>
% <img src="../docs/pub_heated_brick_2.jpg">
% </td></tr>
% <tr><td>Figure 1. Definition of the geometry of the heated brick</td></tr>
% </table>
% </html>

%%
%
% Analytical solution is due to the Reference
% M. Necati Ozisik “Boundary Value Problems of Heat Conduction”.  Dover Publications,
% INC., N.Y. 1989.

%% Solution
%
% The solution flow is very similar to that of pub_heated_brick.  The only
% difference is due to the use of quadratic tetrahedra.  The
% changes in the code will be clearly marked with _[[TETRAHEDRA]]_.
function pub_cfin_tet
    pu= physical_units_struct;
    [fens,fes]=Ansys_mesh_import('cfin.cdb',pu.MM);
    
    %     gv=drawmesh({fens,mesh_boundary(cat(cat(fes{3},fes{1}),fes{2}))},'fes','facecolor','g','shrink',0.75)
    
    %%
    % Define the material properties.
    % Spreader (copper)
    KXX1=3.8600000E-001*pu.W/pu.MM/pu.K;
    % CPU (silicon)
    KXX2=2.2750000E-004*pu.W/pu.MM/pu.K;
    % Fins (aluminum)
    KXX3=1.6700000E-001*pu.W/pu.MM/pu.K;
    h= 3.21658E-006*pu.W/pu.MM^2/pu.K;
    Q=0.0025*pu.W/pu.MM^3;
    %%
    % We are going to create the objects for the analysis.  The thermal
    % property and the thermal material objects:
    prop1=property_heat_diffusion(struct('thermal_conductivity',KXX1));
    mater1=material_heat_diffusion (struct('property',prop1));
    femm1 = femm_heat_diffusion (struct ('material',mater1,...
        'fes',fes{1},...
        'integration_rule',tet_rule(struct( 'npts',4))));
    prop2=property_heat_diffusion(struct('thermal_conductivity',KXX2));
    mater2=material_heat_diffusion (struct('property',prop2));
    femm2 = femm_heat_diffusion (struct ('material',mater2,...
        'fes',fes{2},...
        'integration_rule',tet_rule(struct( 'npts',4))));
    prop3=property_heat_diffusion(struct('thermal_conductivity',KXX3));
    mater3=material_heat_diffusion (struct('property',prop3));
    femm3 = femm_heat_diffusion (struct ('material',mater3,...
        'fes',fes{3},...
        'integration_rule',tet_rule(struct( 'npts',4))));
    
    %%
    % The geometry nodal field is created from the finite element node set.
    geom = nodal_field(struct('name',['geom'], 'dim', 3, 'fens',fens));
    %%
    % The temperature field has one degree of freedom per node.
    temp=nodal_field(struct('name',['temp'], 'dim', 1, 'nfens',geom.nfens));
    
    %%
    % The essential boundary conditions are applied next.
    %%
    % Essential boundary condition: zero temperature on face ABCD.
    %     fenids=fenode_select(fens,struct('box',[-inf,inf,-inf,inf,-36.33,-36.33],...
    %         'inflate', 1/10)) ;
    %     temp = set_ebc(temp, fenids, true, [], 100.0);
    temp = apply_ebc (temp);
    
    
    %%
    % Number the free degrees of freedom
    temp = numberdofs (temp);
    
    amb = 0*temp;
    amb.fixed_values = amb.fixed_values + 20;
    amb = apply_ebc (amb);
    
    %%
    % Calculate and assemble the conductivity matrix.
    K = conductivity(femm1, sysmat_assembler_sparse, geom, temp)...
        + conductivity(femm2, sysmat_assembler_sparse, geom, temp)...
        + conductivity(femm3, sysmat_assembler_sparse, geom, temp);
    
    bfes=mesh_boundary(cat(cat(fes{3},fes{1}),fes{2}));
    bfemm = femm_heat_diffusion (struct ('material',mater3,...
        'fes',subset(bfes,fe_select(fens,bfes,struct('box',[-inf,inf,-inf,inf,-25.4*pu.MM,inf]))),'surface_transfer',h,...
        'integration_rule',tri_rule(struct( 'npts',3))));
    %     drawmesh({fens,bfemm.fes},'fes','facecolor','g')
    K = K + surface_transfer(bfemm, sysmat_assembler_sparse, geom, temp);
    %%
    % Calculate and assemble the non-zero-EBC heat load.
    fi= force_intensity(struct('magn',Q));
    F = surface_transfer_loads(bfemm, sysvec_assembler, geom, temp, amb)...
        + distrib_loads(femm2, sysvec_assembler, geom, temp, fi, 3);
    %%
    % Compute the values of the temperature for the free degrees of freedom and
    % distribute the results to the nodal field.
    temp = scatter_sysvec(temp, K\F);
    [min(temp.values),max(temp.values)]
    %%
    % Graphical rendering  of the computed temperatures. First we show the
    % temperature displayed with filled surfaces on the boundary of the mesh
    % color-coded  with the temperature.
    %%
    % Create the graphics viewer and reset the view. Note: so that we don't have to
    % see all the intermediate figures in the published output, we create
    % the figure as invisible, and then we turn it on when we call interact().
    close all
    figure('visible', 'off')
    gv=graphic_viewer;
    gv=reset (gv,struct('axes',gca));
    
    %%
    % Create the data colormap object for mapping temperatures to colors.
    dcm=data_colormap(struct('range',[min(temp.values),max(temp.values)],'colormap',cadcolors));
    %%
    % Create the color field, with colors at nodes.  The colors will be
    % interpolated using the finite element basis functions on the mesh.
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
        map_data(dcm, temp.values)));
    
    %%
    % Plot the color field on the surface of the brick.
    draw(mesh_boundary(femm1.fes, []), gv, struct ('x',geom, 'u',0*geom,...
        'colorfield',colorfield, 'shrink',1));
    draw(mesh_boundary(femm2.fes, []), gv, struct ('x',geom, 'u',0*geom,...
        'colorfield',colorfield, 'shrink',1));
    draw(mesh_boundary(femm3.fes, []), gv, struct ('x',geom, 'u',0*geom,...
        'colorfield',colorfield, 'shrink',1));
    
    %%
    % Add the color bar
    draw_colorbar(gv, struct('colormap',dcm.colormap,'label','Temperature',...
        'minmax',[min(temp.values),max(temp.values)+eps]));
    %%
    %  Add the labels to the axes...
    labels ([])
    %%
    % ...and make the figure visible and interact with the view.
    interact(gv);
    
    
    %     %%
    %     % We are going to present the temperatures using isosurfaces.  These are
    %     % surfaces of constant temperature, which are going to be color-coded using
    %     % a map from temperatures to colors.
    %     %%
    %     % Create the graphics viewer and reset the view.
    %     close all
    %     figure('visible', 'off')
    %     gv=graphic_viewer;
    %     gv=reset (gv,struct('axes',gca));
    %
    %     %%
    %     % Create the data colormap object for mapping temperatures to colors.
    %     dcm=data_colormap(struct('range',[min(temp.values),max(temp.values)],'colormap',hot));
    %     %
    %     %     %%
    %     %     % Plot the boundary surface of the brick in wireframe rendering for reference.
    %     %     draw(mesh_boundary (bfes, []), gv, struct ('x',geom, 'u',0*geom,...
    %     %         'facecolor','none'));
    %     %     for isovalue = [0.5:0.25:3.0]
    %     %         draw_isosurface(fes,  gv, struct ('x',geom, 'u',0*geom,...
    %     %             'scalarfield',temp,'isovalue',isovalue,'color',map_data(dcm, isovalue)));
    %     %     end
    %     %     %%
    %     %     % Add the color bar
    %     %     draw_colorbar(gv, struct('colormap',dcm.colormap,'label','Temperature',...
    %     %         'minmax',[min(temp.values),max(temp.values)+eps]));
    %     %     %%
    %     %     %  Add the labels to the axes...
    %     %     labels ([])
    %     %
    %     %     %%
    %     %     % Add lighting to aid in the interpretation of the scene.
    %     %     headlight(gv);
    %     %     %%
    %     %     % ...and make the figure visible and interact with the view.
    %     %     interact(gv);
    %     %
    %
    %
    %     %%
    %     % Temperature at point P: Find the node at that point using a bounding box,...
    %     Pid=fenode_select(fens,struct('box',bounding_box(P),...
    %         'inflate', 1/1000)) ;
    %
    %     %%
    %     % ...  and retrieve the temperature:
    %     temp.values(Pid)
    %     %%
    %     % This may be compared with the analytical solution of 0.901 degrees at
    %     % point P. With the given mesh we are therefore within a percent of the
    %     % analytical solution.
    %
    %
    %
    %     %%
    %     % Temperature at point P: Find the node at that point using a bounding box,...
    %     Pid=fenode_select(fens,struct('box',bounding_box(P),...
    %         'inflate', 1/1000)) ;
    %
    %     %%
    %     % ...  and retrieve the temperature:
    %     temp.values(Pid)
    %     %%
    %     % This may be compared with the analytical solution of 0.901 degrees at
    %     % point P. With the given mesh we are therefore within a  percent of the
    %     % analytical solution.
    %
    %
    %
    %
    %     %% Discussion
    %     %
    %     %%
    %     % The quadratic tetrahedron is seen to be quite accurate. If you'd
    %     % like to see the effect of the approximation order on the tetrahedron
    %     % on the quality of the solution, switch to linear tetrahedra (T4) by
    %     % removing (commenting out) the conversion
    %     %
    %     %     [fens,fes]=T4_to_T10(fens,fes);
    %     %
    % %%
    % % For reasonable accuracy one would then have to refine the mesh, and not
    % % only once but twice, to reduce the error below 10%.
end