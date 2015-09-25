% Axially symmetric model of temperature distribution in a concrete column
% submerged in ice water.
alltim = timetic(); tic(alltim);
kappa=1.8*[1.0 0; 0 1.0]; % conductivity matrix
Q = 4.5;
a = 2.5;
t = 1.0;
boundaryf =@(x)0;
N=10000;
tim = timetic(); tic(tim);
[fens,fes]=Q4_block(a,t,N,1, struct('axisymm',true));
     edge_fes = mesh_boundary (fes,struct('axisymm',true));
     toc(tim)
     
tim = timetic(); tic(tim);
% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.fes= fes;
region.reduction='axisymm';
region.Q= Q;
region.integration_rule = gauss_rule(struct('order',2,'dim', 2));
model_data.region{1} =region;

clear essential
essential.temperature=@(x)0;
Outer=fe_select(fens,edge_fes,struct('facing' ,true, 'direction', [1,0]))  ;
essential.fes = subset(edge_fes, Outer);
model_data.boundary_conditions.essential{1} = essential;
 
% Solve
model_data =heat_diffusion_steady_state(model_data);
toc(tim)

max(gather_sysvec(model_data.temp))
%
% model_data.postprocessing.z_scale = 1;
% % model_data.postprocessing.camera =[ 325.7121 -164.7175   21.8267   24.2132   12.8789   -0.6138         0         0    1.0000    7.3068];
% model_data =heat_diffusion_plot_raised_surface(model_data);
%
toc(alltim)