% Layered-square simulation. Varied orthotropic and isotropic layers.

kappaLayer1=[2.25 0; 0 0.06]; % orthotropic conductivity matrix
kappaLayer2=[0.25 0; 0 0.25]; % isotropic conductivity matrix
alpha =45;% local material orientation angle
ca=cos(2*pi/360*alpha); sa=sin(2*pi/360*alpha);
Rm = [ca, -sa;sa, ca];% local material directions
num_integ_pts=1; % 1-point quadrature

% Generate the mesh
[fens,fes, groups] = targe2_mesher({...
    ['curve 1 line -40 20 40 20'],...
    ['curve 2 line 40 20 40 40'],...
    ['curve 3 line 40 40 -40 40'],...
    ['curve 4 line -40 40 -40 20'],...
    ['curve 5 line -40 20 -40 0'],...
    ['curve 6 line -40 0 40 0'],...
    ['curve 7 line 40 0 40 20'],...
    ['curve 8 line -40 0 -40 -20'],...
    ['curve 9 line -40 -20 40 -20'],...
    ['curve 10 line 40 -20 40 0'],...
    ['curve 11 line -40 -20 -40 -40'],...
    ['curve 12 line -40 -40 40 -40'],...
    ['curve 13 line 40 -40 40 -20'],...
    ['subregion 1  property 1 ' ...
    '   boundary 1 2 3 4'],...
    ['subregion 2  property 2 '...
    '   boundary 6 7 -1 5 '],...
    ['subregion 3  property 1 ' ...
    '   boundary 9 10 -6 8'],...
    ['subregion 4  property 2 ' ...
    '   boundary 12 13 -9 11'],...
    ['m-ctl-point constant 2.75']
    }, 1.0);


% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappaLayer1;
region.fes= subset(fes,groups{1});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =Rm;
model_data.region{1} =region;

clear region
region.conductivity =kappaLayer2;
region.fes= subset(fes,groups{2});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =[];
model_data.region{2} =region;

clear region
region.conductivity =kappaLayer1;
region.fes= subset(fes,groups{3});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =Rm;
model_data.region{3} =region;

clear region
region.conductivity =kappaLayer2;
region.fes= subset(fes,groups{4});
region.integration_rule =tri_rule(struct('npts',1));
region.Rm =[];
model_data.region{4} =region;

clear essential
essential.temperature=20;
essential.node_list = [fenode_select(fens,struct('box',[-40 -40 -40 40],...
    'inflate', 0.01))];
model_data.boundary_conditions.essential{1} = essential;

clear essential
essential.temperature=80;
essential.node_list = [fenode_select(fens,struct('box',[40 40 -40 40],...
    'inflate', 0.01))];
model_data.boundary_conditions.essential{2} = essential;

% Solve
model_data =heat_diffusion_steady_state(model_data);


model_data.postprocessing.z_scale = 1;
heat_diffusion_plot_raised_surface(model_data);
set(gca,'zlim', [0,80]);
