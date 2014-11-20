% T4NAFEMS benchmark problem  solved with quadratic elements.
kappa=[52 0; 0 52]; % conductivity matrix
Q=0.0; % uniform heat source
h=750;

mesh_size  =  0.1;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    ['curve 1 line 0 0 0.6 0'],...
    ['curve 2 line 0.6 0 0.6 0.2'],...
    ['curve 3 line 0.6 0.2 0.6 1.0'],...
    ['curve 4 line 0.6 1.0 0 1.0'],...
    ['curve 5 line 0 1.0 0 0'],...
    'subregion 1  property 1 boundary 1 2 3 4 5',...
    ['m-ctl-point constant '  num2str(mesh_size)]
    }, 1.0, struct('quadratic',true));
tolerance =mesh_size/100;

% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.fes= fes;
region.integration_rule =tri_rule(struct('npts',3));
model_data.region{1} =region;

clear convection
convection.ambient_temperature=0;
convection.surface_transfer_coefficient  =h;
convection.fes = subset(edge_fes,...
    [fe_select(fens,edge_fes,struct('box',[0.6 0.6 0 1],'inflate', tolerance)),...
    fe_select(fens,edge_fes,struct('box',[0 1 1 1],'inflate', tolerance))]);
convection.integration_rule=gauss_rule(struct('dim',1,'order',3));
model_data.boundary_conditions.convection{1} = convection;

clear essential
essential.temperature=100;
essential.fes = subset(edge_fes,...
    fe_select(fens,edge_fes,struct('box',[0. 0.6 0 0],'inflate', tolerance)));
model_data.boundary_conditions.essential{1} = essential;


% Solve
model_data =heat_diffusion_steady_state(model_data);

% Report Computed  temperature
disp( ['Temperature at point A '  ])
gather_values(model_data.temp,fenode_select(fens,...
    struct('box',[0.6 0.6 0.2 0.2],'inflate', tolerance)))

model_data.postprocessing.z_scale = 0.01;
heat_diffusion_plot_raised_surface(model_data);
