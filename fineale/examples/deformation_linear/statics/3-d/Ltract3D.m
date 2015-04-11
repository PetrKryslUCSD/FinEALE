% disp('L-shaped domain (one quarter of a square domain with a square hole)');
 
function Ltract_energy
disp('L-shaped domain, 3-D geometry');
% Parameters:
E=1e7;
nu=0.24999999;
graphics=0;
scale = 1;
n=2;

% Mesh
[fens,fes] = Q4_L2x2;
for n=1:1:n
    [fens,fes]=Q4_refine(fens,fes);
end
% [fens,fes]=Q4_to_Q16(fens,fes, []); 
fens.xyz=[0*fens.xyz(:,1),fens.xyz]; % make into a 3-D geometry




% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.E =E;
region.nu =nu;
region.reduction ='stress';
region.fes= fes;
region.integration_rule = gauss_rule (struct('dim', 2, 'order', 2));
region.Rm =[[0, 1, 0]', [0, 0, 1]'];
model_data.region{1} =region;


% first out of plane
clear essential
essential.component= 1;
essential.fixed_value= 0;
essential.node_list =(1:count(fens));;
model_data.boundary_conditions.essential{1} = essential;

% now in plane
clear essential
essential.component= 2;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 0 0 0 0 1],'inflate', 0.0001));
model_data.boundary_conditions.essential{2} = essential;

clear essential
essential.component= 3;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 0 0 1 0 0],'inflate', 0.0001));
model_data.boundary_conditions.essential{3} = essential;

clear essential
essential.component= 3;
essential.fixed_value= 0.25;
essential.node_list = fenode_select (fens,struct('box',[0 0 0 1 1 1],'inflate', 0.0001));
model_data.boundary_conditions.essential{4} = essential;

% Solve
model_data =deformation_linear_statics(model_data);

model_data.postprocessing.u_scale= scale;
model_data.postprocessing.boundary_only= 0;
% model_data=deformation_plot_deformation(model_data);
model_data.postprocessing.output='Cauchy'; 
model_data.postprocessing.component=3; 
model_data.postprocessing.outputRm=eye(3); 
model_data=deformation_plot_stress(model_data);