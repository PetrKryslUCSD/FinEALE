% Illustration of boundary conditions.
kappa=[5.5 0; 0 0.2]; % conductivity matrix
magn = -1;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    ['curve 1 line 0 0 48 0'],...
    ['curve 2 line 48 0 48 48'],...
    ['curve 3 line 48 48 0 48'],...
    ['curve 4 line 0 48 0 26'],...
    ['curve 5 line 0 26 0 13'],...
    ['curve 6 line 0 13 0 0'],...
    'subregion 1  property 1 boundary 1 2 3 4 5 6',...
    ['m-ctl-point constant 1.5']
    }, 1.0);
     
% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.fes= fes;
region.integration_rule = tri_rule(struct('npts',1));
model_data.region{1} =region;

clear essential
essential.temperature=0;
essential.fes = subset(edge_fes,[edge_groups{[1, 2, 3, 4, 6]}]);
model_data.boundary_conditions.essential{1} = essential;
 
clear flux
flux.normal_flux=magn;
flux.fes = subset(edge_fes,[edge_groups{5}]);
flux.integration_rule =simpson_1_3_rule(struct('dim',1));
model_data.boundary_conditions.flux{1} = flux;

% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 4;
model_data.postprocessing.camera =[ 325.7121 -164.7175   21.8267   24.2132   12.8789   -0.6138         0         0    1.0000    7.3068];
mmodel_data =heat_diffusion_plot_raised_surface(model_data);
