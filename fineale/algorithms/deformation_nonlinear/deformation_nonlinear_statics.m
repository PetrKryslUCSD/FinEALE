function model_data = deformation_nonlinear_statics(model_data)
% Algorithm for static nonlinear deformation (stress) analysis.
%
% function model_data = deformation_nonlinear_statics(model_data)
%
% Arguments
% model_data = struct  with fields as follows.
%
% model_data.fens = finite element node set (mandatory)
%
% For each region (connected piece of the domain made of a particular material),
% mandatory:
% model_data.region= cell array of struct with the attributes, each region
%           gets a struct with attributes
%     fes= finite element set that covers the region
%     integration_rule =integration rule
%     material = material object,
%
%        If the material orientation matrix is not the identity and needs
%        to be supplied,  include the attribute
%     Rm= constant orientation matrix or a handle  to a function to compute
%           the orientation matrix (see the class femm_base).
%
% Example:
%     clear region
%     region.material=mater;
%     region.fes= fes;
%     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
%     model_data.region{1} =region;
%
% Example:
%     clear region
%     region.material=mater;
%     region.fes= fes;% set of finite elements for the interior of the domain
%     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
%     region.Rm =@LayerRm;
%     model_data.region{1} =region;
%
% Instead of the  finite element  set, the material, and the integration
% rule, one can supply directly the FEMM.
%     femm =  the finite element model machine
%
% Example:
%     clear region
%     prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
%     region.femm= femm_deformation_nonlinear_h8msgso(...
%         struct ('material',material_deformation_stvk_triax(struct('property',prop)),...
%         'fes',fes, ...
%         'integration_rule',gauss_rule(struct('dim',3,'order',2))));;
%     model_data.region{1} =region;
%
% For essential boundary conditions (optional):
% model_data.boundary_conditions.essential = cell array of struct,
%           each piece of surface with essential boundary condition gets one
%           element of the array with a struct with the attributes
%           as recognized by the set_ebc() method of nodal_field
%     is_fixed=is the degree of freedom being prescribed or
%           freed? (boolean)
%     component=which component is affected (if not supplied,
%           or if supplied as empty default is all components)
%     fixed_value=fixed (prescribed) displacement (scalar, array, or function handle)
%           function v=fixedValue(Lambda)
%           Lambda=load parameter value,
%           v =two-dimensional array, one row for each node on the list,
%           one column for each component specified.
%
%     fes = finite element set on the boundary to which
%                       the condition applies
%               or alternatively
%     node_list = list of nodes on the boundary to which
%                       the condition applies
%           Only one of fes and node_list needs to be given.
%
% Example:
%     clear essential
%     essential.component= [1,3];
%     essential.fixed_value= 0;
%     essential.node_list = [[fenode_select(fens, struct('box', [0,a,0,0,-Inf,Inf],...
%         'inflate',tolerance))],[fenode_select(fens, struct('box', [0,a,b,b,-Inf,Inf],...
%         'inflate',tolerance))]];;
%     model_data.boundary_conditions.essential{2} = essential;
%
% Example:
%     clear essential
%     essential.component= [1];
%     essential.fixed_value= 0;
%     essential.fes = mesh_boundary(fes);
%     model_data.boundary_conditions.essential{1} = essential;
%
% For traction boundary conditions (optional):
% model_data.boundary_conditions.traction = cell array of struct,
%           each piece of surface with traction boundary condition gets one
%           element of the array with a struct with the attributes
%     traction=traction (vector), supply a zero for component in which
%           the boundary condition is inactive
%     fes = finite element set on the boundary to which
%                       the condition applies
%     integration_rule= integration rule
%
% Example:
%     clear traction
%     traction.fes =subset(bdry_fes,bclx);
%     traction.traction= @(x) ([sigmaxx(x);sigmaxy(x);0]);
%     traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));%<==
%     model_data.boundary_conditions.traction{1} = traction;
%
% For body loads (optional):
% model_data.body_load = cell array of struct,
%          each piece of the domain can have each its own body load
%     force  = force density vector
%     fes = finite element set to which the load applies
%     integration_rule= integration rule
%
%
% Control parameters:
% The following attributes  may be supplied as fields of the model_data struct:
%     load_multipliers = For what load multipliers should dissolution be
%           calculated? Array of monotonically increasing numbers.
%     renumber = true or false flag (default is true)
%     renumbering_method = optionally choose the renumbering
%           method  ('symrcm' or 'symamd')
%     line_search = Should we use line search? Boolean.  Default = true.
%     maxdu_tol = Tolerance on the magnitude  of the largest incremental 
%           displacement component.
%     iteration_observer = handle to an observer function to be called 
%           after each iteration was computed.  Default is  to do nothing.
%           The observer function has a signature
%                     function iteration_observer(lambda,iter,du,model_data)
%           where lambda is the current load factor, iter is the iteration 
%           number, du is the nodal field of current displacement increments.
%     increment_observer = handle of an observer function 
%           to be called after convergence is reached in each step (optional)
%           The observer function has a signature
%                     function output(lambda, model_data)
%           where lambda is the current load factor
%
% Output
% model_data =  structure that was given on input updated with:
% The struct model_data on output incorporates in addition to the input the fields
%     geom = geometry field
%     u = computed displacement nodal field
%     reactions = computed reaction nodal field
%


renumber = true; % Should we renumber?
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

line_search = ~true; % Should we use line search?
if (isfield(model_data,'line_search'))
    line_search  =model_data.line_search;
end

load_multipliers = 1.0; % For what load multipliers should dissolution be calculated?
if (isfield(model_data,'load_multipliers'))
    load_multipliers  =model_data.load_multipliers;
end

iteration_observer =[];%  do nothing
if ( isfield(model_data,'iteration_observer'))
    iteration_observer  =model_data.iteration_observer;;
end

increment_observer =@(lambda,model_data) disp(['Load multiplier = ' num2str(lambda)]);
if ( isfield(model_data,'load_increment_observer'))% This is an alias
    increment_observer  =model_data.load_increment_observer;;
end
if ( isfield(model_data,'increment_observer'))
    increment_observer  =model_data.increment_observer;;
end

maxdu_tol=0;% Tolerance on the magnitude  of the largest incremental displacement component
if (isfield(model_data,'maxdu_tol'))
    maxdu_tol  =model_data.maxdu_tol;
end

maxbal_tol=0;% Tolerance on the magnitude  of the out-of-balance force
if (isfield(model_data,'maxbal_tol'))
    maxbal_tol  =model_data.maxbal_tol;
end

% Extract the nodes
fens =model_data.fens;

% Construct the geometry field
geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));
model_data.geom = geom;
            
% Construct the displacement field
un1 = nodal_field(struct('name',['un1'], 'dim', geom.dim, 'nfens', geom.nfens));

% Apply the essential boundary conditions on the displacement field
if (isfield(model_data.boundary_conditions, 'essential' ))
    for j=1:length(model_data.boundary_conditions.essential)
        essential =model_data.boundary_conditions.essential{j};
        if (isfield( essential, 'fes' ))
            essential.node_list= connected_nodes(essential.fes);
        end
        is_fixed=ones(length(essential.node_list),1);
        if (isfield( essential, 'is_fixed' ))
            is_fixed= essential.is_fixed;
        end
        essential.is_fixed=is_fixed;
        component=[];% Implies all components
        if (isfield( essential, 'component' ))
            component= essential.component;
        end
        if (isempty(component))
            component =1:un1.dim;
        end
        fixed_value =0;
        % If the essential boundary condition is not  supplied by a
        % function, set it equal to zero initially.
        if (isfield( essential, 'fixed_value' )) && (~isa(essential.fixed_value,'function_handle'))
            fixed_value = essential.fixed_value;
        end
        val=zeros(length(essential.node_list),1)+fixed_value;
        for k=1:length( component)
            un1 = set_ebc(un1, essential.node_list, is_fixed, component(k), val);
        end
        un1 = apply_ebc (un1);
        % If the essential boundary condition is to be load-factor-dependent, the
        % values must be supplied by a function.   Find out.
        essential.load_factor_dependent =(isa(essential.fixed_value,'function_handle'));
        model_data.boundary_conditions.essential{j} =essential;
    end
    clear essential fenids fixed component fixed_value  val
end

% Number the equations
un1 = numberdofs (un1, Renumbering_options);
un =un1;% Converged displacement at the beginning of the current step, i. e.  u_{n}
unm1 =un;% Converged displacement one step back, i. e.  u_{n-1}

% Create the necessary FEMMs
for i=1:length(model_data.region)
    region =model_data.region{i};
    % Create  the finite element model machine, if not supplied
    if (~isfield( region, 'femm'))
        mater = region.material;% Material needs to be defined
        Rm = [];%  Is the material orientation matrix supplied?
        if (isfield( region, 'Rm'))
            Rm= region.Rm;
        end
        region.femm = femm_deformation_nonlinear (struct ('material',mater,...
            'fes',region.fes,...
            'integration_rule',region.integration_rule,'Rm',Rm));
    end
    % Give the  FEMM a chance  to precompute  geometry-related quantities
    region.femm = associate_geometry(region.femm,geom);
    model_data.region{i}=region;
    clear region prop mater Rm  femm
end

% Initially the stiffness matrix is empty.  We can check it  in the
% treatment of nonzero essential boundary conditions and constructed if
% needed.
K =[];

% Increment magnitudes
load_increments=diff(unique(sort([0,load_multipliers])));

% For all  load increments
for incr =1:length(load_increments) % Load-implementation loop
    
    lambda =sum(load_increments(1:incr));
    dlambda=load_increments(incr);
    
    % Initial value of the displacement  in the current step
    un1 =un;
    
    % Process load-factor-dependent essential boundary conditions:
    aany_nonzero_EBC = false;
    if ( isfield(model_data,'boundary_conditions')) && ...
            (isfield(model_data.boundary_conditions, 'essential' ))
        for j=1:length(model_data.boundary_conditions.essential)
            if (model_data.boundary_conditions.essential{j}.load_factor_dependent)% this needs to be computed
                essential =model_data.boundary_conditions.essential{j};
                fixed_value =essential.fixed_value(lambda);
                for k=1:length(essential.component)
                    un1 = set_ebc(un1, essential.node_list, ...
                        essential.is_fixed, essential.component(k), fixed_value(:,k));
                    aany_nonzero_EBC = aany_nonzero_EBC || any(fixed_value(:,k)~=0);
                end
                un1 = apply_ebc (un1);
            end
            clear essential
        end
    end
    
    % Initialization
    un1 = apply_ebc(un1); % Apply EBC
    du = 0*un1; % Displacement increment
    du = apply_ebc(du);
    un1Allfree=un1;% This field will hold the displacements for all notes and degrees of freedom.
    un1Allfree.is_fixed(:)=0;% all displacements will be free
    un1Allfree = numberdofs (un1Allfree);
        
    % Initialize the load vector
    F =zeros(un1.nfreedofs,1);
    
    % If any boundary conditions are inhomogeneous, calculate  the force
    % vector due to the displacement increment. Then update the guess of
    % the new displacement.
    if (aany_nonzero_EBC)
        for i=1:length(model_data.region)
            F = F + nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, geom, un, unm1, un1-un, dlambda);
        end
        % Provided we got converged results  in the last step, we already
        % have a usable stiffness matrix.
        if (isempty(K)) % we don't have a stiffness matrix
            K=  sparse(un1.nfreedofs,un1.nfreedofs);
            for i=1:length(model_data.region)
                K = K + stiffness(model_data.region{i}.femm, sysmat_assembler_sparse, geom, un, unm1, dlambda);
                K = K + stiffness_geo(model_data.region{i}.femm, sysmat_assembler_sparse, geom, un, unm1, dlambda);
            end
        end
        un1 = un1 + scatter_sysvec(du, K\F);
    end
    
        
    % Iteration loop
    iter=1;
    while true %  Iteration loop
        
        % Initialize the load vector
        F =0*F;
        
        % Construct the system stiffness matrix. 
        K=  sparse(un1.nfreedofs,un1.nfreedofs);
        for i=1:length(model_data.region)
            K = K + stiffness(model_data.region{i}.femm, sysmat_assembler_sparse, geom, un1, un, dlambda);
            K = K + stiffness_geo(model_data.region{i}.femm, sysmat_assembler_sparse, geom, un1, un, dlambda);
        end
        
        % Process the body load
        if (isfield(model_data, 'body_load' ))
            for j=1:length(model_data.body_load)
                body_load =model_data.body_load{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',body_load.fes,...
                    'integration_rule',body_load.integration_rule));
                fi= force_intensity(struct('magn',body_load.force));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 3);
            end
            clear body_load fi  femm
        end
        
        % Process the traction boundary condition
        if (isfield(model_data.boundary_conditions, 'traction' ))
            for j=1:length(model_data.boundary_conditions.traction)
                traction =model_data.boundary_conditions.traction{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',traction.fes,...
                    'integration_rule',traction.integration_rule));
                fi= force_intensity(struct('magn',traction.traction));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 2);
            end
            clear traction fi  femm
        end
        
        % Process the nodal force boundary condition
        if (isfield(model_data.boundary_conditions, 'nodal_force' ))
            for j=1:length(model_data.boundary_conditions.nodal_force)
                nodal_force =model_data.boundary_conditions.nodal_force{j};
                femm = femm_deformation (struct ('material',[],...
                    'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
                    'integration_rule',point_rule));
                fi= force_intensity(struct('magn',nodal_force.force));
                F = F + distrib_loads(femm, sysvec_assembler, geom, un1, fi, 0);
            end
            clear nodal_force fi femm
        end
        
        %  External loads vector
        FL = lambda * F;
        
        % The restoring  force vector
        FR =zeros(un1.nfreedofs,1);
        for i=1:length(model_data.region)
            FR = FR + restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1, un, dlambda); 
        end

        % Solve the system of linear algebraic equations
        F = (FL + FR);
        dusol= K\F;
        
        %  Distribute the solution
        du = scatter_sysvec(du, dusol);
        
        % Do we need line search?
        eta= 1.0;
        if (line_search)
            R0 = dot(F,dusol);
            for  linesrch   =1:1 % How many  searches should we do?
                FR =zeros(un1.nfreedofs,1);
                for i=1:length(model_data.region)
                    FR = FR + restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1+eta*du, un, dlambda);
                end
                F = FL + FR;
                R1 = dot(F,dusol);
                a = R0/R1;
                if ( a<0 )
                    eta = a/2 +sqrt((a/2)^2 -a);
                else
                    eta =a/2;
                end
                eta=min( [eta, 1.0] );
            end
        end
        
        % Increment the displacements
        un1 = un1 + eta*du;  
        
        % Compute the increment data
        maxdu  =eta*max(max(abs(du.values)));
        maxbal=max(abs(F));
        
        % Report  iteration results?
        if (~isempty(iteration_observer))
            model_data.un = un;
            model_data.un1 = un1;
            model_data.dt = dlambda;
            model_data.maxdu = maxdu;
            model_data.maxbal = maxbal;
            iteration_observer (lambda,iter,du,model_data);
        end
        
        % Converged?
        if (maxdu <= maxdu_tol)  || (maxbal <= maxbal_tol)
            break; % Breakout of the iteration loop
        end;    
        
        iter=iter+1;
        
    end %  Iteration loop
    
     
    % Update the model data
    model_data.un = un;
    model_data.un1 = un1;
    model_data.dt = dlambda;
    %  Now the reactions
    un1Allfree.values=un1.values;% This field holds the current converged displacements;
    unAllfree=un1Allfree;
    unAllfree.values=un.values;% This field holds the current converged displacements;
    % since all the degrees of freedom are free, the restoring forces
    % can be used to compute the reactions. Note the negative sign: the
    % reactions are the opposite of the resisting forces of the
    % material.
    FR =zeros(un1Allfree.nfreedofs,1);
    for i=1:length(model_data.region)% Compute the resisting forces of the material
        FR = FR - restoring_force(model_data.region{i}.femm, sysvec_assembler, geom, un1Allfree, unAllfree, dlambda);
    end
    model_data.reactions=scatter_sysvec(un1Allfree,FR);
        
    % Converged
    % Update the FEMMs
    for i=1:length(model_data.region)
       [~,model_data.region{i}.femm]  =...
           restoring_force(model_data.region{i}.femm, sysvec_assembler, ...
           geom, un1, un, dlambda); 
    end
    
    % Report results
   if ~isempty(increment_observer)% report the progress
        increment_observer (lambda,model_data);
    end
    
    % Reset  the displacement field for the next load step
    unm1 =un;
    un = un1;
    F =0*F;
    
end % Load-incrementation loop

end
