% # poisson_fem.jl
% #
% # Main file for solving the 2D Poisson equation over the unit square.
% # Depends on finite_element.jl
% #
% # Amuthan Arunkumar Ramabathiran
% # (aparyap@gmail.com)
% #
% # Distributed under The Code Project Open License (CPOL)
% # http://www.codeproject.com/info/cpol10.aspx
alltim = timetic(); tic(alltim);
kappa=[1.0 0; 0 1.0]; % conductivity matrix
magn = -6;
boundaryf =@(x)1 + x(1)^2 + 2*x(2)^2;
N=10;
tim = timetic(); tic(tim);
[fens,fes]=T3_block(1.0,1.0,N,N, 1.0);
     edge_fes = mesh_boundary (fes,  []);
     toc(tim)
     
tim = timetic(); tic(tim);
% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.fes= fes;
region.Q= magn;
region.integration_rule = tri_rule(struct('npts',1));
model_data.region{1} =region;

clear essential
essential.temperature=@(x)(1 + x(:,1).^2 + 2*x(:,2).^2);
essential.fes = edge_fes;
model_data.boundary_conditions.essential{1} = essential;
 
% Solve
model_data =heat_diffusion_steady_state(model_data);
toc(tim)

%
model_data.postprocessing.z_scale = 1;
% model_data.postprocessing.camera =[ 325.7121 -164.7175   21.8267   24.2132   12.8789   -0.6138         0         0    1.0000    7.3068];
model_data =heat_diffusion_plot_raised_surface(model_data);

toc(alltim)