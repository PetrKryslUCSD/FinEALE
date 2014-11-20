% L-shaped region with uniform heat source.
% L-shaped region (one quarter of a rectangular region with a rectangular
% hole).   Uniform heat source, temperature fixed on the boundary.
%
kappa=[0.2 0; 0 0.2]; % conductivity matrix
Q=0.013; % uniform heat source
[fens,fes] = targe2_mesher({...
    ['curve 1 line 20 0 48 0'],...
    ['curve 2 line 48 0 48 48'],...
    ['curve 3 line 48 48 0 48'],...
    ['curve 4 line 0 48 0 13'],...
    ['curve 5 line 0 13 20 13'],...
    ['curve 6 line 20 13 20 0'],...
    'subregion 1  property 1 boundary 1 2 3 4 5 6',...
    ['m-ctl-point constant 2.5']
    }, 1.0);

% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.Q =Q;
region.fes= fes;
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =[];
model_data.region{1} =region;

clear essential
essential.temperature=30;
essential.node_list = [fenode_select(fens,struct('box',[48 48 0 48],...
    'inflate', 0.01)),fenode_select(fens,struct('box',[0 48 48 48],...
    'inflate', 0.01)),fenode_select(fens,struct('box',[0 20 0 13],...
    'inflate', 0.01))];
model_data.boundary_conditions.essential{1} = essential;

model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 1;
model_data.postprocessing.camera =[-214.1067  -64.2844  201.6982   19.8741   23.1974   13.4603    0.5637    0.2108    0.7986 10.3396];
model_data.postprocessing.draw_mesh = true;
model_data =heat_diffusion_plot_raised_surface(model_data);

model_data.postprocessing.scale = 35;
model_data=heat_diffusion_plot_integration_points(model_data);
set(gca,'zlim', [0,50]);