function t4nafems_Q8_mech
% T4NAFEMS benchmark problem  solved with quadratic elements.
pu= physical_units_machine;
kappa=[52 0; 0 52]*pu('W/kg/m'); % conductivity matrix
Q=0.0*pu('W/m^3'); % uniform heat source
h=750*pu('W/kg/m^2');
E=206.8*pu('GPa');;
nu= 0.3;
CTE= 7.5e-6;
nel=4;

mesh_size  =  0.03*pu('m');
gtol=mesh_size/100;
[fens,fes]=Q8_block((0.6)*pu('m'),(1.0)*pu('m'),2*nel,5*nel,1.0*pu('m'));

% Thermal solution
% =========================================================================
% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.conductivity =kappa;
region.fes= fes;
region.integration_rule =gauss_rule(struct('dim',2,'order',2));
model_data.region{1} =region;

clear convection
convection.ambient_temperature=0;
convection.surface_transfer_coefficient  =h;
edge_fes= mesh_boundary(fes, struct('other_dimension', 1.0*pu('m')));
ix=[fe_select(fens,edge_fes,struct('box',[0.6 0.6 0 1]*pu('m'),'inflate', gtol)),...
    fe_select(fens,edge_fes,struct('box',[0. 0.6 1 1]*pu('m'),'inflate', gtol))];
convection.fes = subset(edge_fes,ix);
convection.integration_rule=gauss_rule(struct('dim',1,'order',3));
model_data.boundary_conditions.convection{1} = convection;

clear essential
essential.temperature=100;
essential.fes = subset(edge_fes,...
    fe_select(fens,edge_fes,struct('box',[0. 0.6 0 0],'inflate', gtol)));
model_data.boundary_conditions.essential{1} = essential;


% Solve
model_data =heat_diffusion_steady_state(model_data);

% Report Computed  temperature
disp( ['Temperature at point A '  ])
gather_values(model_data.temp,fenode_select(fens,...
    struct('box',[0.6 0.6 0.2 0.2],'inflate', gtol)))

model_data.postprocessing.z_scale = 0.01;
heat_diffusion_plot_raised_surface(model_data);

% Temperature field
Temp=model_data.temp;
% Temp.values(:)=100;% Set the temperature increment to a constant value

pause(2)

% Mechanical solution
% =========================================================================
% Compose the model data
clear model_data
model_data.fens =fens;

clear region
region.E =E;
region.nu =nu;
region.alpha =CTE;
region.reduction ='strain';
region.fes= fes;
region.integration_rule =gauss_rule(struct('dim',2,'order',2));
region.Rm =[];
model_data.region{1} =region;

clear essential % the symmetry plane
essential.component= 1;
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[0 0 -inf inf],'inflate',gtol));
model_data.boundary_conditions.essential{1} = essential;


clear essential % the symmetry plane
essential.component= [1,2];
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct('box',[-inf inf 0 0],'inflate',gtol));
model_data.boundary_conditions.essential{2} = essential;

% Thermal loading through this temperature field
model_data.temperature=Temp;

% Solve
model_data =deformation_linear_statics(model_data);

disp([num2str(max( model_data.u.values ))])

model_data.postprocessing.u_scale=260;
model_data.postprocessing.quantity='U2';
model_data.postprocessing.boundary_only=false;
model_data.postprocessing.cmap=jet(12);
model_data=deformation_plot_deformation(model_data);
view (2)
end