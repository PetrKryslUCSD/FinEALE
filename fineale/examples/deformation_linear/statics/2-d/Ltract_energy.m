% disp('L-shaped domain (one quarter of a square domain with a square hole)');
 
function Ltract_energy
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


% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.E =E;
region.nu =nu;
region.reduction ='strain';
region.fes= fes;
region.integration_rule = gauss_rule (struct('dim', 2, 'order', 2));
region.Rm =[];
model_data.region{1} =region;

clear essential
essential.component= 1;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 0 0 1],'inflate', 0.0001));
model_data.boundary_conditions.essential{1} = essential;

clear essential
essential.component= 2;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 1 0 0],'inflate', 0.0001));
model_data.boundary_conditions.essential{2} = essential;

clear essential
essential.component= 2;
essential.fixed_value= 0.25;
essential.node_list = fenode_select (fens,struct('box',[0 1 1 1],'inflate', 0.0001));
model_data.boundary_conditions.essential{3} = essential;

% Solve
model_data =deformation_linear_statics(model_data);

model_data.postprocessing.u_scale= scale;
model_data.postprocessing.boundary_only= 0;
model_data=deformation_plot_deformation(model_data);