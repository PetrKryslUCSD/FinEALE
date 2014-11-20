% T4NAFEMS benchmark problem  solved with linear and quadratic elements.
%  Convergence study.
kappa=[52 0; 0 52]; % conductivity matrix
Q=0.0; % uniform heat source
h=750;
online_graphics =true;
num_integ_pts=[1,6];
mesh_options = [struct('quadratic',false),struct('quadratic',true)];
results = {[], []};
% mesh_sizes = [1.0, 0.5, 0.25, 0.125 0.125/2 0.125/4 0.125/8  0.125/16];
mesh_sizes = [1.0, 0.5, 0.25, 0.125 0.125/2 0.125/4];
% mesh_sizes = [1.0, 0.5, 0.25, 0.125 0.125/2];
for et= 1:2
    for mesh_size = mesh_sizes
        
        [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
            ['curve 1 line 0 0 0.6 0'],...
            ['curve 2 line 0.6 0 0.6 0.2'],...
            ['curve 3 line 0.6 0.2 0.6 1.0'],...
            ['curve 4 line 0.6 1.0 0 1.0'],...
            ['curve 5 line 0 1.0 0 0'],...
            'subregion 1  property 1 boundary 1 2 3 4 5',...
            ['m-ctl-point constant ' num2str(mesh_size)],...
            ['m-ctl-point 1 xy 0.6 0.2 near ' num2str(mesh_size/10)...
            ' influence ' num2str(mesh_size/16)],...
            ['m-ctl-point 1 xy 0.6 0 near ' num2str(mesh_size/10)...
            ' influence ' num2str(mesh_size/8)]
            }, 1.0, mesh_options(et));
        tolerance =mesh_size/100;
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.conductivity =kappa;
        region.fes= fes;
        region.integration_rule =tri_rule(struct('npts',num_integ_pts(et)));
        model_data.region{1} =region;
        
        clear convection
        convection.ambient_temperature=0;
        convection.surface_transfer_coefficient  =h;
        convection.fes = subset(edge_fes,...
            [fe_select(fens,edge_fes,struct('box',[0.6 0.6 0 1],'inflate', tolerance)),...
            fe_select(fens,edge_fes,struct('box',[0 1 1 1],'inflate', tolerance))]);
        convection.integration_rule=gauss_rule(struct('dim',1,'order',4));
        model_data.boundary_conditions.convection{1} = convection;
        
        clear essential
        essential.temperature=100;
        essential.fes = subset(edge_fes,...
            fe_select(fens,edge_fes,struct('box',[0. 0.6 0 0],'inflate', tolerance)));
        model_data.boundary_conditions.essential{1} = essential;
        
        
        % Solve
        model_data =heat_diffusion_steady_state(model_data);
        
        results{et} =[results{et},gather_values(model_data.temp,fenode_select(fens,...
            struct('box',[0.6 0.6 0.2 0.2],'inflate', tolerance)))];
        
        % Plot
        if online_graphics
            model_data.postprocessing.z_scale = 0.01;
            heat_diffusion_plot_raised_surface(model_data);
            pause (1)
        end
    end
    results{et}
end
xs =results{2};% Estimate from the quadratic triangle results
[xestim, beta] = richextrapol(xs(3:end),mesh_sizes(3:end));
set(0, 'DefaultAxesFontSize', 18);
figure
loglog(mesh_sizes,abs(results{2}-xestim)/xestim,'bo-','linewidth',3)
hold on
grid on
loglog(mesh_sizes,abs(results{1}-xestim)/xestim,'rs-','linewidth',3)
xlabel('log(mesh size)')
ylabel('log(|estimated temperature error|)')
set_graphics_defaults

figure
loglog(mesh_sizes(2:end),abs(diff(results{2})),'bo-','linewidth',3)
hold on
loglog(mesh_sizes(2:end),abs(diff(results{1})),'rs-','linewidth',3)
grid on
xlabel('log(mesh size)')
ylabel('log(|approximate temperature error|)')
set_graphics_defaults
