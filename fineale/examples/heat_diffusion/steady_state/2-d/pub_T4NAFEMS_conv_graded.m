%% Two-dimensional heat transfer with convection:  study with graded meshes
%

%%
% Link to the  <matlab:edit('pub_T4NAFEMS_conv_graded') m-file>.
%

%% Description
%
% Consider a plate of uniform thickness, measuring 0.6 m by 1.0 m. On one
% short edge the temperature is fixed at 100 °C, and on one long edge the
% plate is perfectly insulated so that the heat flux is zero through that
% edge. The other two edges are losing heat via convection to an ambient
% temperature of 0 °C. The thermal conductivity of the plate is 52.0 W/(m
% .°K), and the convective heat transfer coefficient is 750 W/(m^2.°K).
% There is no internal generation of heat. Calculate the temperature 0.2 m
% along the un-insulated long side, measured from the intersection with the
% fixed temperature side. The reference result is 18.25 °C.

%%
%
% <html> <table border=0><tr><td> <img src="../docs/pub_T4NAFEMS_conv.jpg"
% width="70%"> </td></tr> <tr><td>Figure 1. Definition of the geometry of
% the domain</td></tr> </table> </html>


%%
% The reference temperature at the point A  is 18.25 °C according to the
% NAFEMS publication ( hich cites the book Carslaw, H.S. and J.C. Jaeger,
% Conduction of Heat in Solids. 1959: Oxford University Press).

%%
% The present  tutorial will investigate the reference temperature  and it
% will attempt to  estimate the  limit value more precisely using a
% sequence of graded meshes and Richardson's extrapolation.

%% Solution
%
function  pub_T4NAFEMS_conv_graded
    pu=physical_units_struct;%  Bring in  definitions of physical units
    kappa=[52 0; 0 52]*pu.W/(pu.M*pu.K); % conductivity matrix
    h=750*pu.W/(pu.M^2*pu.K);% surface heat transfer coefficient
    Width=0.6*pu.M;% Geometrical dimensions
    Height=1.0*pu.M;
    HeightA=0.2*pu.M;
    Thickness=1.0*pu.M;
    Refinement_factors = 1./2.^(-1:1:4);% Refinement factors for the sequence of simulations
    tolerance =HeightA*min(Refinement_factors)/1000;
    
    
    %%
    % The simulation will be executed inside the loop over all the mesh
    % sizes. The array  |results| will collect the temperatures at point A.
    results = [];
    for Refinement_factor = Refinement_factors
        
        %%
        % Generate the triangle mesh for the current mesh size. The domain
        % is given as a sequence of vertices starting at the lower left
        % corner, traversing the boundary in counterclockwise sense.  Edge
        % 1 is between the first and the second vertex and so on. This call
        % creates a cell array of commands that define the curves of the
        % boundary.
        Commands =  targe2_mesher_vl_cmds([0,0; Width,0;  Width,HeightA;  ...
            Width,Height;   0,Height; ]);
        %%
        % Define the "subregion" by specifying its boundary curves (this
        % corresponds to region  in our model_data).
        Commands{end+1}='subregion 1  property 1 boundary 1 2 3 4 5';
        
        %%
        % Now specify the mesh size. The "background" mesh size, which is
        % the mesh size desired everywhere where the mesh-size controls  to
        % be defined next do not apply.
        mesh_size  = Refinement_factor*Height;
        Commands{end+1}=['m-ctl-point constant ' num2str(mesh_size)];
        
        %%
        % We define two mesh-size control points.  The first at the
        % location  of the point A. The desired mesh size is 1/10 of the
        % background mesh size,  and the influence region is twice the
        % desired size of the elements
        local_mesh_size  =Refinement_factor*Height/10;
        Commands{end+1}=['m-ctl-point 1 xy ' num2str([Width,HeightA]) ...
            ' near ' num2str(local_mesh_size) ...
            ' influence ' num2str(local_mesh_size*2)];
        
        %%
        % The second mesh-size control point is located at the point of the
        % singularity at the transition between the prescribed temperature
        % and the convection boundary condition.
        local_mesh_size  =Refinement_factor*Height/10;
        Commands{end+1}=['m-ctl-point 1 xy ' num2str([Width,0])...
            ' near ' num2str(local_mesh_size)...
            ' influence ' num2str(local_mesh_size*2)];
        
        %%
        % Note that the mesh size refinement described by the variable
        % |Refinement_factor| varies inside the loop.   However, the
        % proportions of the background mesh size and the desired mesh
        % sizes  at the target point and the singularity are fixed.   This
        % is important for the Richardson's extrapolation which relies  on
        % there being a single refinement factor across the entire mesh:
        % the ratio of the mesh size  in the mesh (j+1) to the mesh size in
        % the mesh (j) at any point within the domain  must be a constant.
        
        %%
        % The automatic mesher can now be run to triangulate the domain
        % using the Commands information. Note that we are requesting
        % quadratic triangles to be generated.
        [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher(...
            Commands,Thickness,struct('quadratic',true));
        
        %%
        % Set up the model data.  The nodes:
        clear model_data
        model_data.fens =fens;
        
        %%
        % The region: note our use of six point quadrature  for the
        % quadratic triangles.
        clear region
        region.conductivity =kappa;
        region.fes= fes;
        region.integration_rule =tri_rule(struct('npts',6));
        model_data.region{1} =region;
        
        
        %%
        % The convection boundary condition is applied along the edges
        % 2,3,4. The elements along the boundary are quadratic line
        % elements L3. The order-four Gauss quadrature is sufficiently
        % accurate.
        clear convection
        convection.ambient_temperature=0;
        convection.surface_transfer_coefficient  =h;
        convection.fes = subset(edge_fes,[edge_groups{2},edge_groups{3},edge_groups{4}]);
        convection.integration_rule=gauss_rule(struct('dim',1,'order',4));
        model_data.boundary_conditions.convection{1} = convection;
        
        %%
        % The prescribed temperature is applied along edge 1 (the bottom
        % edge in Figure 1)..
        clear essential
        essential.temperature=100;
        essential.fes = subset(edge_fes,[edge_groups{1}]);
        model_data.boundary_conditions.essential{1} = essential;
        
        
        %%
        % The model data is defined, solve for the temperatures.
        model_data =heat_diffusion_steady_state(model_data);
        
        %%
        % Collect the temperature  at the point A  [coordinates
        % (Width,HeightA)].
        results =[results,gather_values(model_data.temp,fenode_select(fens,...
            struct('box',[  Width,Width,HeightA,HeightA],'inflate', tolerance)))];
        
        
        %%
        % Plot the temperature as a raised surface.
        model_data.postprocessing.z_scale = 0.01;
        heat_diffusion_plot_raised_surface(model_data);
        
    end
    
    %%
    % These are the computed results for the temperature at point A:
    results
    
    %%
    % Richardson extrapolation is used to estimate the true solution from
    % the results for the finest three meshes.
    [xestim, beta] = richextrapol(results(end-2:end),Refinement_factors(end-2:end));
    disp(['Estimated true solution for temperature at A: ' num2str(xestim) ' degrees'])
    
    %%
    % Plot the estimated true error.
    figure
    loglog(Refinement_factors,abs(results-xestim)/xestim,'bo-','linewidth',3)
    grid on
     xlabel('log(refinement factor)')
    ylabel('log(|estimated temperature error|)')
    set_graphics_defaults
    
    %%
    % The estimated true error has  a slope of approximately 3 on the
    % log-log scale. This is more closely aligned with the theoretical
    % prediction of the convergence rate for quadratic elements.
    %%
    % Plot the absolute values of the approximate error (differences  of
    % successive solutions).
    figure
    loglog(Refinement_factors(2:end),abs(diff(results)),'bo-','linewidth',3)
    grid on
    xlabel('log(refinement factor)')
    ylabel('log(|approximate temperature error|)')
    set_graphics_defaults
    
    
%% Discussion
% 
    %%
    % The use of graded mesh-size meshes is  more efficient than use of
    % uniform meshes. The extrapolation seems to be well justified, and
    % provides a good quality estimate of the true solution.
end