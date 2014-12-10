function model_data =heat_diffusion_steady_state(model_data)
% Steady-state heat conduction solver.
%
% function model_data =heat_diffusion_steady_state(model_data)
%
%
% Arguments: 
%
% model_data = struct  with fields as follows.
%
% model_data.fens = finite element node set (mandatory)
%
% For each region (connected piece of the domain made of a particular material),
% mandatory:
% model_data.region= cell array of struct with the attributes, each region 
%           gets a struct with attributes
%     region.conductivity = material conductivity
%     region.Q = material internal heat generation rate
%     region.fes= finite element set that covers the region
%     region.integration_rule =integration rule
%     
% For essential boundary conditions (optional):
% model_data.boundary_conditions.essential = cell array of struct,
%           each piece of surface with essential boundary condition gets one
%           element of the array with a struct with the attributes
%     essential.temperature=fixed (prescribed) temperature (scalar),  or
%           handle to a function with signature
%               function T =f(x)
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
% For flux boundary conditions (optional):
% model_data.boundary_conditions.flux = cell array of struct,
%           each piece of surface with flux boundary condition gets one
%           element of the array with a struct with the attributes
%     flux.normal_flux=ambient temperature (scalar)
%     flux.fes = finite element set on the boundary to which 
%                       the condition applies
%     flux.integration_rule= integration rule
%    
% Control parameters:
% model_data.renumber = true or false flag (default is true)
% model_data.renumbering_method = optionally choose the renumbering 
%       method  (symrcm or symamd)
%
% Output:
% model_data= the struct on input is augmented with 
%     for each region and for each boundary condition of the convection or 
%         flux type the finite element model is added to model_data.  For 
%         instance model_data.region{1}.fem or, if convection BC is defined, 
%         model_data.boundary_conditions.convection{1}.fem 
%     renumbering information, if renumbering was requested      
%     geom =the nodal field that is the geometry
%     temp =the nodal field that is the computed temperature 
    
%     Control parameters
    renumber = ~true; % Should we renumber?
    if (isfield(model_data,'renumber'))
        renumber  =model_data.renumber;
    end
    Renumbering_options =struct( [] );
    
    % Should we renumber the nodes to minimize the cost of the solution of 
    % the coupled linear algebraic equations?
    if (renumber)
        renumbering_method = 'symamd'; % default choice
        if ( isfield(model_data,'renumbering_method'))
            renumbering_method  =model_data.renumbering_method;;
        end
        % Run the renumbering algorithm
        model_data =renumber_mesh(model_data, renumbering_method);;
        % Save  the renumbering  (permutation of the nodes)
        clear  Renumbering_options; Renumbering_options.node_perm  =model_data.node_perm;
    end
    
    % Extract the nodes
    fens =model_data.fens;
    
    % Construct the geometry field
    geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));
    
    % Construct the temperature field
    temp = nodal_field(struct('name',['temp'], 'dim', 1, 'nfens', geom.nfens));
    
    % Apply the essential boundary conditions on the temperature field
    if (isfield(model_data.boundary_conditions, 'essential' ))
        for j=1:length(model_data.boundary_conditions.essential)
            essential =model_data.boundary_conditions.essential{j};
            if (isfield( essential, 'fes' ))
                fenids= connected_nodes(essential.fes);
            else
                fenids= essential.node_list;
            end
            fixed=ones(length(fenids),1);
            comp=[];
            T_fixed=0;
            if (isfield( essential, 'temperature' ))
                if (strcmp(class(essential.temperature),'function_handle'))
                   T_fixed = essential.temperature(geom.values(fenids,:));
                else
                    T_fixed = essential.temperature;
                end
            end
            if (length(T_fixed)==1), T_fixed =repmat(T_fixed,length(fenids),1); end
            val=zeros(length(fenids),1)+T_fixed;
            temp = set_ebc(temp, fenids, fixed, comp, val);
            temp = apply_ebc (temp);
        end
        clear essential fenids fixed comp T_fixed  val
    end
    
    % Number the equations
    temp = numberdofs (temp,Renumbering_options);
    
    % Initialize the heat loads vector
    F =zeros(temp.nfreedofs,1);
    
    % Create the heat diffusion finite element models for the regions
    for i=1:length(model_data.region)
        region =model_data.region{i};
        % Construct the property  and material  objects
        prop=property_heat_diffusion(struct('thermal_conductivity',region.conductivity));
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
    if (isfield(model_data.boundary_conditions, 'convection' ))
        amb = 0*clone(temp,'amb'); % create the ambient temperature field
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
    if (isfield(model_data.boundary_conditions, 'flux' ))
        for j=1:length(model_data.boundary_conditions.flux)
            flux =model_data.boundary_conditions.flux{j};
            flux.femm = femm_heat_diffusion (struct ('material',[],...
                'fes',flux.fes,...
                'integration_rule',flux.integration_rule));
            model_data.boundary_conditions.flux{j}=flux;    
        end
        clear flux fi  femm
    end
     
     
    % Construct the system conductivity matrix 
    K=  sparse(temp.nfreedofs,temp.nfreedofs); % (all zeros, for the moment)
    for i=1:length(model_data.region)
        region =model_data.region{i};
        % Add up all the conductivity matrices for all the regions        
        K = K + conductivity(region.femm, sysmat_assembler_sparse, geom, temp); 
        Q = 0;% default is no internal heat generation rate
        if (isfield( region, 'Q'))%  Was it defined?
            Q= region.Q;
        end
        if (Q~=0)% If it was supplied, and it is nonzero, compute its contribution.
            fi= force_intensity(struct('magn',Q));
            F = F + distrib_loads(region.femm, sysvec_assembler, geom, temp, fi, 3);
        end
        % Loads due to the essential boundary conditions on the temperature field
        if (isfield(model_data.boundary_conditions, 'essential' ))
            F = F + nz_ebc_loads_conductivity(region.femm, sysvec_assembler, geom, temp);
        end
        clear region Q prop mater Rm 
    end
     
    % Process the convection boundary condition
    if (isfield(model_data.boundary_conditions, 'convection' ))
        amb = 0*clone(temp,'amb'); % create the ambient temperature field
        for j=1:length(model_data.boundary_conditions.convection)
            convection =model_data.boundary_conditions.convection{j};
            % Apply the fixed ambient temperature
            fenids= connected_nodes(convection.fes);
            fixed=ones(length(fenids),1);
            comp=[];
            T_fixed =0;% default
            if (isfield( convection, 'ambient_temperature' ))
                T_fixed = convection.ambient_temperature;
            end
            val=zeros(length(fenids),1)+T_fixed;
            % This looks like we are setting essential boundary conditions, 
            % but in reality we are modifying the ambient temperature
            amb = set_ebc(amb, fenids, fixed, comp, val);
            amb = apply_ebc (amb);
            K = K + surface_transfer(convection.fem, sysmat_assembler_sparse, geom, temp);   
            F = F + surface_transfer_loads(convection.fem, sysvec_assembler, geom, temp, amb); 
            % Note that EBC will contribute through the surface heat transfer matrix            
            if (isfield(model_data.boundary_conditions, 'essential' ))
                F = F + nz_ebc_loads_surface_transfer(convection.fem, sysvec_assembler, geom, temp);
            end
        end
        clear convection fenids fixed comp T_fixed  val 
    end
    
    % Process the flux boundary condition
    if (isfield(model_data.boundary_conditions, 'flux' ))
        for j=1:length(model_data.boundary_conditions.flux)
            flux =model_data.boundary_conditions.flux{j};
            fi= force_intensity(struct('magn',flux.normal_flux));
            % Note the sign  which reflects the formula (negative sign
            % in front of the integral)
            F = F - distrib_loads(flux.femm, sysvec_assembler, geom, temp, fi, 2);   
        end
        clear flux fi  
    end
     
    % Solve for the temperatures
    temp = scatter_sysvec(temp, K\F);
    
    % Update the model data
    model_data.geom = geom;
    model_data.temp = temp;
end
