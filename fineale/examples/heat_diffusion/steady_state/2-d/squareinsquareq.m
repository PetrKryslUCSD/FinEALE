% Square-in-square, heat flux visualization.
kappainner=[2.25 0; 0 0.06]; % orthotropic conductivity matrix
kappaouter=[0.25 0; 0 0.25]; % isotropic conductivity matrix
alpha =-45;% local material orientation angle
ca=cos(2*pi/360*alpha); sa=sin(2*pi/360*alpha);
Rm = [ca, -sa;sa, ca];% local material directions

[fens,fes, groups] = targe2_mesher({...
    ['curve 1 line -48 -48 48 -48'],...
    ['curve 2 line 48 -48 48 48'],...
    ['curve 3 line 48 48 -48 48'],...
    ['curve 4 line -48 48 -48 -48'],...
    ['curve 5 line 0 -31 31 0'],...
    ['curve 6 line 31 0 0 31'],...
    ['curve 7 line 0 31 -31 0'],...
    ['curve 8 line -31 0 0 -31'],...
    ['subregion 1  property 1 ' ...
    '   boundary 1 2 3 4 -8 -7 -6 -5'],...
    ['subregion 2  property 2 '...
    '   boundary 5 6 7 8'],...
    ['m-ctl-point constant 4.75']
    }, 1.0);

% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappainner;
region.fes= subset(fes,groups{2});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =Rm;
model_data.region{1} =region;

clear region
region.conductivity =kappaouter;
region.fes= subset(fes,groups{1});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =[];
model_data.region{2} =region;

clear essential
essential.temperature=20;
essential.node_list = [fenode_select(fens,struct('box',[-48 48 -48 -48],...
    'inflate', 0.01))];
model_data.boundary_conditions.essential{1} = essential;

clear essential
essential.temperature=57;
essential.node_list = [fenode_select(fens,struct('box',[-48 48 48 48],...
    'inflate', 0.01))];
model_data.boundary_conditions.essential{2} = essential;

% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.scale = 45;
model_data=heat_diffusion_plot_integration_points(model_data);
model_data.postprocessing.z_scale = 1;
mmodel_data=heat_diffusion_plot_level_curves(model_data);
% heat_diffusion_plot_raised_surface(model_data, options);
set(gca,'zlim', [0,60]);
view (2);


