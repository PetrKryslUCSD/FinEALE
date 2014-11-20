% Analogy between heat conduction and shear stress distribution in an U-shape 

kappa=[0.2 0; 0 0.2]; % conductivity matrix
Q=0.01; % uniform heat source
num_integ_pts=1; % 1-point quadrature
W=150; H=260; tw=90; th=120;% some dimensions
%  Mesh the domain
[fens,fes] = targe2_mesher_vl([...
    -W,0; W,0; W,H;  W-tw,H; W-tw,th; -W+tw,th; -W+tw,H; -W,H; ...
    ], 1.0,struct('mesh_size', tw/6 ));
 

% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.Q =Q;
region.fes= fes;
region.integration_rule = tri_rule(struct('npts',1));
model_data.region{1} =region;

clear essential
essential.temperature=0;
essential.fes = mesh_boundary(fes,[]);
model_data.boundary_conditions.essential{1} = essential;
 
% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 1;
figure;
model_data =heat_diffusion_plot_raised_surface(model_data);
model_data.postprocessing.gv=[];

figure; 
model_data =heat_diffusion_plot_level_curves(model_data);
model_data.postprocessing.gv=[];

model_data.postprocessing.color = 'k';
model_data.postprocessing.draw_mesh= false;;
figure; 
model_data =heat_diffusion_plot_raised_surface(model_data);
model_data =heat_diffusion_plot_level_curves(model_data);
