% Chimney with convection boundary condition.
function ChimneyN_conv
kappa=0.5*eye(2); % conductivity matrix
hi=5;
ho=15;
Q=0; % uniform heat source
a=1;
Ti = 1000;
To = -20;
mesh_sizes= 0.5./2.^(0:1:4);
maxT=[];
for mesh_size=mesh_sizes
    % Generate the mesh of the domain
    [fens,fes,groups,edge_fes,edge_groups] =  targe2_mesher({...
        ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(2*a) ' ' num2str(0) ''],...
        ['curve 2 line ' num2str(2*a) ' ' num2str(0) ' ' num2str(a) ' ' num2str(a) ''],...
        ['curve 3 line ' num2str(a) ' ' num2str(a) ' ' num2str(0) ' ' num2str(a) ''],...
        ['curve 4 line ' num2str(0) ' ' num2str(a) ' ' num2str(0) ' ' num2str(0) ''],...
        'subregion 1  property 1 boundary 1 2 3 4',...
        ['m-ctl-point constant ' num2str( mesh_size )]
        }, 1.0);
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.conductivity =kappa;
    region.Q = 0;
    region.fes= fes;
    region.integration_rule =tri_rule(struct('npts', 1));
    model_data.region{1} =region;
    
    clear convection
    convection.ambient_temperature=Ti;
    convection.surface_transfer_coefficient  =hi;
    convection.fes =subset(edge_fes,[edge_groups{3}]);
    convection.integration_rule=gauss_rule(struct('dim',1,'order',2));
    model_data.boundary_conditions.convection{1} = convection;
    
    clear convection
    convection.ambient_temperature=To;
    convection.surface_transfer_coefficient  =ho;
    convection.fes =subset(edge_fes,[edge_groups{1}]);
    convection.integration_rule=gauss_rule(struct('dim',1,'order',2));
    model_data.boundary_conditions.convection{2} = convection;
    
    % Solve
    model_data =heat_diffusion_steady_state(model_data);
    
    maxT= [maxT,max(model_data.temp.values)];
end
maxT
diff(maxT)
A= [log(mesh_sizes(end-2:end))',ones(3,1)];
b=log(abs(diff(maxT(end-3:end))))';
A\b
assignin('caller','fineale_test_passed',((norm(9.048450848159243e2-maxT(end)'))<1e-2));