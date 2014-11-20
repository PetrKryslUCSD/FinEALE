function model_data =heat_diffusion_transient(model_data)
% Transient linear heat conduction solver. 
%
% function model_data =heat_diffusion_transient(model_data)
%
% The generalized 
% trapezoidal rule is used for direct time integration.
%
% Input: model is defined by the struct as follows.
%
% model_data = struct with attributes
% model_data.fens = finite element node set (mandatory)
%
% For each region (connected piece of the domain made of a particular material),
% mandatory:
% model_data.region= cell array of struct with the attributes, each region 
%           gets a struct with attributes
%     region.conductivity = material conductivity
%     region.specific_heat = material specific heat per unit volume
%     region.fes= finite element set that covers the region
%     region.integration_rule =integration rule
%     
% For essential boundary conditions (optional):
% model_data.boundary_conditions.essential = cell array of struct,
%           each piece of surface with essential boundary condition gets one
%           element of the array with a struct with the attributes
%     essential.temperature=fixed (prescribed) temperature (scalar)
%     essential.fes = finite element set on the boundary to which 
%                       the condition applies
%               or alternatively
%     essential.node_list = list of nodes on the boundary to which 
%                       the condition applies
%               Only one of essential.fes and essential.node_list needs a given.
%            
% For convection boundary conditions (optional):
% model_data.boundary_conditions.convection = cell array of struct,
%           each piece of surface with convection boundary condition gets one
%           element of the array with a struct with the attributes
%     convection.ambient_temperature=ambient temperature (scalar)
%     convection.surface_transfer_coefficient  = surface heat transfer coefficient
%     convection.fes = finite element set on the boundary to which 
%                       the condition applies
%     convection.integration_rule= integration rule
%    
% For initial condition:
% model_data.initial_condition = struct with the attribute
%     initial_condition.temperature=initial temperature as a scalar value or 
%            as a function of the location array geom.values to return the 
%            initial temperature at each node
%        
% Control parameters:
% model_data.theta= parameter of the generalized trapezoidal rule, theta=0 is 
%           the explicit Euler, theta=1 is the implicit Euler, and theta=1/2 is 
%           the Crank-Nicholson, (optional, default 1.0)
% model_data.dt= time step (mandatory)
% model_data.tend= end of the integration interval (mandatory)
% model_data.observer = handle of an observer function (optional).  
%           The observer function is called before the integration starts 
%           and then after the completion of each step.
%           Signature: function observer (t,  model_data)
% model_data.renumber = true or false flag (default is true)
% model_data.renumbering_method = optionally choose the renumbering 
%       method  (symrcm or symamd)
%
%
% Output:
% model_data= the struct on input is augmented with 
%     renumbering information, if renumbering was requested      
%     geom =the nodal field that is the geometry
%     temp =the nodal field that is the computed temperature 
    
%     Control parameters
    theta = 1; % parameter of the generalized trapezoidal method
    if ( isfield(model_data,'theta'))
        theta  =model_data.theta;;
    end
    observer =@(t,model_data) disp(['Time ' num2str(t)]);
    if ( isfield(model_data,'observer'))
        observer  =model_data.observer;;
    end
    dt =model_data.dt; % time step
    tend =model_data.tend; % integrate up to the time tend
    
    renumber = true; % should we renumber the equations?
    if (isfield(model_data,'renumber'))
        renumber  =model_data.renumber;
    end
    Renumbering_options =struct( [] );
    
    % Should we renumber the nodes to minimize the cost of the solution of 
    % the coupled linear algebraic equations?
    if (renumber)
        renumbering_method = 'symamd'; % default choice (for Choleski)
        if ( isfield(model_data,'renumbering_method'))
            renumbering_method  =model_data.renumbering_method;;
        end
        % Run the renumbering algorithm
        model_data =renumber_mesh(model_data, renumbering_method);;
        % Save  the renumbering  (permutation of the nodes)
        clear  Renumbering_options; Renumbering_options.node_perm  =model_data.node_perm;
    end
     
%     Extract the nodes
    fens =model_data.fens;
    
%     Construct the geometry field
    geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));
    
%     Construct the temperature field
    tempn = nodal_field(struct('name',['tempn'], 'dim', 1, 'nfens', geom.nfens));
    
%     Apply the essential boundary conditions on the temperature field
    if ( isfield(model_data,'boundary_conditions')) && (isfield(model_data.boundary_conditions, 'essential' ))
        for j=1:length(model_data.boundary_conditions.essential)
            essential =model_data.boundary_conditions.essential{j};
            if (isfield( essential, 'fes' ))
                essential.fenids= connected_nodes(essential.fes);
            else
                essential.fenids= essential.node_list;
            end
            essential.fixed=ones(length(essential.fenids),1);
            essential.comp=[];
            val=zeros(length(essential.fenids),1)+0;
            tempn = set_ebc(tempn, essential.fenids, essential.fixed, essential.comp, val);
            tempn = apply_ebc (tempn);
            model_data.boundary_conditions.essential{j} =essential;
        end
        clear essential fenids fixed comp T_fixed  val
    end
    
%     Number the equations
    tempn = numberdofs (tempn,Renumbering_options);
    
%     Initialize the heat loads vector
    F =zeros(tempn.nfreedofs,1);
    
    % Create the heat diffusion finite element models for the regions
    for i=1:length(model_data.region)
        region =model_data.region{i};
        % Construct the property  and material  objects
        prop=property_heat_diffusion(struct('thermal_conductivity',region.conductivity,...
            'specific_heat',region.specific_heat));
        mater=material_heat_diffusion (struct('property',prop));
        Rm = [];
        if (isfield( region, 'Rm'))
            Rm= region.Rm;
        end
        % This is the model  for the current region: note that we supply 
        % integration  rule and the  material orientation matrix
        region.femm = femm_heat_diffusion (struct ('material',mater,...
                'fes',region.fes,...
                'integration_rule',region.integration_rule,'Rm',Rm));
        model_data.region{i} =region;
        clear region  prop mater Rm 
    end
    
    % Create the finite element models for the boundary conditions
    % Convection boundary conditions:
    if isfield(model_data,'boundary_conditions') && ...
            (isfield(model_data.boundary_conditions, 'convection' ))
        for j=1:length(model_data.boundary_conditions.convection)
            convection =model_data.boundary_conditions.convection{j};
            % Note that we need to supply the surface heat transfer coefficient
            convection.fem = femm_heat_diffusion (struct ('material',[],...
                'fes',convection.fes,...
                'integration_rule',convection.integration_rule,...
                'surface_transfer', convection.surface_transfer_coefficient));
            model_data.boundary_conditions.convection{j} =convection;    
        end
        clear convection 
    end
    %  Flux boundary condition:
    if isfield(model_data,'boundary_conditions') && ...
            (isfield(model_data.boundary_conditions, 'flux' ))
        for j=1:length(model_data.boundary_conditions.flux)
            flux =model_data.boundary_conditions.flux{j};
            flux.fem = femm_heat_diffusion (struct ('material',[],...
                'fes',flux.fes,...
                'integration_rule',flux.integration_rule));
            model_data.boundary_conditions.flux{j}=flux;    
        end
        clear flux   
    end
     
%     Construct the system conductivity and capacity matrix
    K=  sparse(tempn.nfreedofs,tempn.nfreedofs);
    C=  sparse(tempn.nfreedofs,tempn.nfreedofs);
    for i=1:length(model_data.region)
        region =model_data.region{i};
        Q = 0;% default is no internal heat generation rate
        if (isfield( region, 'Q'))%  Was it defined?
            Q= region.Q;
        end
        if (Q~=0)% If it was supplied, and it is nonzero, compute its contribution.
            fi= force_intensity(struct('magn',Q));
            F = F + distrib_loads(region.femm, sysvec_assembler, geom, temp, fi, 3);
        end
        K = K + conductivity(region.femm, sysmat_assembler_sparse, geom, tempn); 
        C = C + capacity(region.femm, sysmat_assembler_sparse, geom, tempn); 
        % %  Loads due to the essential boundary conditions on the temperature field
        % if ( isfield(model_data,'boundary_conditions')) && (isfield(model_data.boundary_conditions, 'essential' ))
            % F = F + nz_ebc_loads_conductivity(femm, sysvec_assembler, geom, tempn);
        % end
    end
     
%     Process the convection boundary condition
    if ( isfield(model_data,'boundary_conditions')) && (isfield(model_data.boundary_conditions, 'convection' ))
        amb = 0*clone(tempn,'amb'); % create the ambient temperature field
        for j=1:length(model_data.boundary_conditions.convection)
            convection =model_data.boundary_conditions.convection{j};
            fenids= connected_nodes(convection.fes);
            fixed=ones(length(fenids),1);
            comp=[];
            T_fixed =0;
            if (isfield( convection, 'ambient_temperature' ))
                T_fixed = convection.ambient_temperature;
            end
            val=zeros(length(fenids),1)+T_fixed;
            amb = set_ebc(amb, fenids, fixed, comp, val);
            amb = apply_ebc (amb);
            K = K + surface_transfer(convection.fem, sysmat_assembler_sparse, geom, tempn);   
            F = F + surface_transfer_loads(convection.fem, sysvec_assembler, geom, tempn, amb);            
            if (isfield(model_data.boundary_conditions, 'essential' ))
                F = F + nz_ebc_loads_surface_transfer(convection.fem, sysvec_assembler, geom, tempn);
            end
        end
        clear convection fenids fixed comp T_fixed  val
    end
    
%     Set the initial condition
    t=0;
    if isfield(model_data,'initial_condition')
        if (strcmp(class(model_data.initial_condition.temperature),'double'))
            vtempn=gather_sysvec(tempn)*0+model_data.initial_condition.temperature;
        else
            tempn.values = model_data.initial_condition.temperature(geom.values);
            vtempn=gather_sysvec(tempn);
        end
    else
        vtempn=gather_sysvec(tempn)*0;    
    end
    tempn = scatter_sysvec(tempn,vtempn);
    model_data.geom = geom;
    
    
    % Factorize the left-hand side matrix for efficiency (Choleski)
    [mA,status,prm] = chol((1/dt*C+theta*K),'lower','vector');%
    if ( status ~= 0 ) error('Choleski factorization failed'), end
    mAt= mA';
    
    %     Solve for the temperatures
    t=0;
    %     First output is the initial condition
    model_data.temp = tempn;% save the computed temperature field
    if ~isempty(observer)% report the progress
        observer (t,model_data);
    end
    
    while true % Time stepping
        t=t+dt;% Advance the current time
        F = 0*F;% time independent part of the temperature load
        Tn=gather_sysvec(tempn);% current temperature system vector
        tempn1 = tempn;% this field will be modified by the boundary conditions
        % Process essential boundary conditions if they are time-dependent
        if ( isfield(model_data,'boundary_conditions')) && ...
                (isfield(model_data.boundary_conditions, 'essential' ))
            for j=1:length(model_data.boundary_conditions.essential)
                essential =model_data.boundary_conditions.essential{j};
                if (isa(essential.temperature,'function_handle'))
                    T_fixed =essential.temperature(t);
                else
                    T_fixed =essential.temperature;
                end
                tempn1 = set_ebc(tempn1, essential.fenids, essential.fixed, essential.comp, T_fixed);
                tempn1 = apply_ebc (tempn1);
            end
            for i=1:length(model_data.region)
                region =model_data.region{i};
                F = F + nz_ebc_loads_conductivity(region.femm, sysvec_assembler, geom, ...
                            (theta*tempn1) + ((1-theta)*tempn));
                F = F + nz_ebc_loads_capacity(region.femm, sysvec_assembler, geom, ...
                            (tempn1-tempn)*(1/dt));
            end
            clear essential  T_fixed   region
        end
        % Needs to be done: process time-dependent convection and time-dependent flux
        % Solve for the new temperature  vector
        bv=((1/dt*C-(1-theta)*K)*Tn+F);% complete right-hand side (loads)
        Tn1(prm) = mAt\(mA\bv(prm)); % compute next temperature
        tempn = scatter_sysvec(tempn1,Tn1);% update the temperature field
        model_data.temp = tempn;% save the computed temperature field
        if ~isempty(observer)% report the progress
            observer (t,model_data);
        end
        if (t==tend)% Have we reached the end?  If so jump out.
            break;
        end
        if (t+dt>tend) % Adjust the last time step so that we exactly reach tend
            dt =tend-t;
        end
    end
     
end




