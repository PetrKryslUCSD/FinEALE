% Orthotropic balloon inflation, axially symmetric model


% Parameters:
E1=1.0;
E2=1.0;
E3=3.0;
nu12=0.29;
nu13=0.29;
nu23=0.19;
G12=0.3;
G13=0.3;
G23=0.3;
p= 0.15;
rin=1;
rex =1.2;
n=2;

% Mesh'
[fens,fes]=Q4_block(rex-rin,pi/2,5,20,struct('axisymm',true));
bdry_fes = mesh_boundary(fes, struct('axisymm', true));
icl = fe_select(fens, bdry_fes, struct('box', [0,0,0,pi/2],'inflate',rin/100));
xy=fens.xyz;
for i=1:count(fens)
    r=rin+xy(i,1); a=xy(i,2);
    xy(i,:)=[r*cos(a) r*sin(a)];
end
fens.xyz=xy;


% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.property = 'orthotropic';
region.E1 =E1;
region.E2 =E2;
region.E3 =E3;
region.G12=G12;
region.G13=G13;
region.G23=G23;
region.nu12=nu12;
region.nu13=nu13;
region.nu23=nu23;
region.reduction ='axisymm';
region.fes= fes;
region.integration_rule = gauss_rule (struct('dim', 2,  'order', 2));
region.Rm =[];
model_data.region{1} =region;

clear essential % the symmetry plane
essential.component= 2;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 rex 0 0],'inflate',rex/10000));
model_data.boundary_conditions.essential{1} = essential;

clear essential % the axis of symmetry
essential.component= 1;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 0 0 rex],'inflate',rex/10000));
model_data.boundary_conditions.essential{2} = essential;

clear traction
traction.fes= subset(bdry_fes,icl);
traction.integration_rule = gauss_rule (struct('dim', 1,  'order', 2));
traction.traction = @(x) (p*x'/norm(x));
model_data.boundary_conditions.traction{1} = traction;


% Solve
model_data =deformation_linear_statics(model_data);

model_data.postprocessing.u_scale= 1;
model_data.postprocessing.stress_component=3;
%model_data.postprocessing.stress_range = 3*[-p,+p];
mmodel_data=deformation_plot_stress(model_data);
view (2);;


