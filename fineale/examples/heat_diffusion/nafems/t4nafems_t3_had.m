% T4NAFEMS benchmark problem adaptively solved with linear elements.

function t4nafems_t3_had
    kappa=[52 0; 0 52]; % conductivity matrix
    Q=0.0; % uniform heat source
    h=750;
    tolerance = 0.2/10000;
    region_definition={'REGEXT 0 0 0.6 1.0',...
        ['curve 1 line 0 0 0.6 0'],...
        ['curve 2 line 0.6 0 0.6 0.2'],...
        ['curve 3 line 0.6 0.2 0.6 1.0'],...
        ['curve 4 line 0.6 1.0 0 1.0'],...
        ['curve 5 line 0 1.0 0 0'],...
        ['subregion 1  property 1 boundary 1 2 3 4 5']};
    
    clear adaptivity_options
    adaptivity_options.mesh_options =struct('quadratic',~true);
    adaptivity_options.targetnel=700;
    adaptivity_options.convergence_rate=1;
    adaptivity_options.initial_mesh_size=0.125;
    adaptivity_options.nadapt=5;
    adaptivity_options.observer =@observer;
    
    function  model_data = make_model_data (fens,fes,groups,edge_fes,edge_groups)
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.conductivity =kappa;
        region.fes= fes;
        region.integration_rule =tri_rule(struct('npts',1));
        model_data.region{1} =region;
        
        clear convection
        convection.ambient_temperature=0;
        convection.surface_transfer_coefficient  =h;
        convection.fes = subset(edge_fes,...
            [fe_select(fens,edge_fes,struct('box',[0.6 0.6 0 1],'inflate', tolerance)),...
            fe_select(fens,edge_fes,struct('box',[0 1 1 1],'inflate', tolerance))]);
        convection.integration_rule=gauss_rule(struct('dim',1,'order',2));
        model_data.boundary_conditions.convection{1} = convection;
        
        clear essential
        essential.temperature=100;
        essential.fes = subset(edge_fes,...
            fe_select(fens,edge_fes,struct('box',[0. 0.6 0 0],'inflate', tolerance)));
        model_data.boundary_conditions.essential{1} = essential;
    end
    
    model_data =heat_diffusion_adaptive_2D_steady_state(region_definition,@make_model_data,adaptivity_options);
    
    % Plot
    function  observer(step,model_data)
        disp(['Adaptive step ' num2str(step), ', ' num2str(count(model_data.region{1}.fes)) '  elements']);
        model_data.postprocessing.z_scale = 0.01;
        heat_diffusion_plot_raised_surface(model_data);
    end
end