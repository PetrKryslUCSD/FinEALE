% Demonstration of under-integration
kappa=[0.2 0; 0 0.2]; % conductivity matrix
Q=0.005; % uniform heat source
radius=48;
n=8;
[fens,fes]=Q4_block(radius,radius,n,n,1.0); 


% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.Q =Q;
region.fes= fes;
region.integration_rule = gauss_rule(struct('dim',2,'order',2));
model_data.region{1} =region;

clear essential
essential.temperature=0;
fenids=connected_nodes(mesh_boundary(fes,[]));
fenids=setdiff(fenids,[fenode_select(fens,struct('box',[0 0.99*radius/2 0 0],...
    'inflate', 0.01*radius/2)),...
    fenode_select(fens,struct('box',[0 0 0 .99*radius/2],...
    'inflate', 0.01*radius/2))]);
essential.node_list = fenids;
model_data.boundary_conditions.essential{1} = essential;


% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 1;
model_data.postprocessing.camera =[ -205.6165 -261.0291  202.6875   13.5376   24.5781   -5.1586 0 0    1. 6.8854];
heat_diffusion_plot_raised_surface(model_data);


% Now change  the  quadrature rule and observe the effect

model_data.region{1}.integration_rule = gauss_rule(struct('dim',2,'order',1));


% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 1;
model_data.postprocessing.camera =[ -205.6165 -261.0291  202.6875   13.5376   24.5781   -5.1586 0 0    1. 6.8854];
heat_diffusion_plot_raised_surface(model_data);
