function model_data = deformation_nonlinear_direct_explicit_CD(model_data)
% Direct-integration explicit solver for mechanical deformation (no damping).
%
% Algorithm for dynamic nonlinear deformation with explicit  (centered-difference)
% direct integration. The structure is assumed to be without damping.
%
% function model_data = deformation_nonlinear_direct_explicit_CD(model_data)
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
%     property = material property hint (optional): 'isotropic' (default), 
%           'orthotropic',...
%        For isotropic property (default):
%     E = material Young's modulus
%     nu = material Poisson's ratio
%        For orthotropic property:
%     E1, E2, E3 = material Young's modulus 
%           in the three material directions
%     G12, G13, G23 = material shear modulus 
%           between the three material directions
%     nu12, nu13, nu23 = material Poisson's ratio
%     rho = mass density (optional for statics, mandatory for dynamics)
%        If the region is for a two-dimensional model (plane strain, plane
%        stress, or axially symmetric) the attribute reduction  needs to be
%        specified.
%     reduction = 'strain' or 'stress' or 'axisymm'
%
%        If the material orientation matrix is not the identity and needs
%        to be supplied,  include the attribute
%     Rm= constant orientation matrix or a handle  to a function to compute
%           the orientation matrix (see the class femm_base).
%
% Example:
%     clear region
%     region.E =E;
%     region.nu=nu;
%     region.fes= fes;
%     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
%     model_data.region{1} =region;
%
% Example:
%     clear region
%     region.property = 'orthotropic';
%     region.rho =rho;
%     region.E1 =E1;     region.E2 =E2;     region.E3 =E3;
%     region.G12=G12;     region.G13=G13;     region.G23=G23;
%     region.nu12=nu12;     region.nu13=nu13;     region.nu23=nu23;
%     region.fes= fes;% set of finite elements for the interior of the domain
%     region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
%     region.Rm =@LayerRm;
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
%     fixed_value=fixed (prescribed) displacement (scalar or function handle)
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
%
% For initial condition:
% model_data.initial_condition = struct with the attribute
%     u_fixed_value=initial displacement as a single
%           vector value (for instance [0,0,-1]) or
%           as a function of the location array geom.values to return the
%           initial displacement at each node; for example
%                 u_fixed_value= f, where
%                     f=@(x)0*x+ones(size(x,1),1)*[-0.02,0,0];
%           Note that x is an array,  one row per node.
%     v_fixed_value=initial velocity as a single
%           vector value (for instance [0,0,-1]) or
%           as a function of the location array geom.values to return the
%           initial velocity at each node; for example
%                 v_fixed_value= f, where
%                     f=@(x)0*x+ones(size(x,1),1)*[0,0,vmag];
%           Note that x is an array,  one row per node.
%
% Example:
%     clear initial_condition
%     initial_condition.u_fixed_value= @(x)0*x;
%     initial_condition.v_fixed_value= @(x)0*x;
%     model_data.initial_condition = initial_condition;
%
% Control parameters:
% The following attributes  may be supplied as fields of the model_data struct:
%      tend= end of the integration interval (mandatory)
%      step_reduction=the algorithm computes the stable time step for
%           explicit integration.  The time step than may be modified (increased)
%           by multiplying it with step_reduction.  Default is step_reduction=1.0.
%      observer = handle of an observer function (optional)
%           The observer function has a signature
%                     function output(t, model_data)
%           where t is the current time
%
% Output
% model_data = updated structure that was given on input updated as:
% It incorporates in addition to the input the nodal_fields
%     geom = geometry field
%     u = displacement field


%     Control parameters
observer =@(t,model_data) disp(['Time ' num2str(t)]);
if ( isfield(model_data,'observer'))
    observer  =model_data.observer;;
end
tend =0;
if ( isfield(model_data,'tend'))
    tend  =model_data.tend;;
end
dt =[];
if ( isfield(model_data,'dt'))
    dt  =model_data.dt;;
end
steps_between_dt_estimations =10;
if ( isfield(model_data,'steps_between_dt_estimations'))
    steps_between_dt_estimations  =model_data.steps_between_dt_estimations;;
end
dt_reduction =10;
if ( isfield(model_data,'dt_reduction'))
    dt_reduction  =model_data.dt_reduction;;
end

% Extract the nodes
fens =model_data.fens;

% Construct the geometry field
model_data.geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));

% Construct the displacement field
model_data.un1 = nodal_field(struct('name',['un1'], 'dim', model_data.geom.dim, 'nfens', model_data.geom.nfens));

% Apply the essential boundary conditions on the displacement field
if isfield(model_data,'boundary_conditions')...
        && (isfield(model_data.boundary_conditions, 'essential' ))
    for j=1:length(model_data.boundary_conditions.essential)
        essential =model_data.boundary_conditions.essential{j};
        if (isfield( essential, 'fes' ))
            essential.fenids= connected_nodes(essential.fes);
        else
            essential.fenids= essential.node_list;
        end
        if (~isfield( essential, 'is_fixed' )) || isempty(essential.is_fixed)
            essential.is_fixed=ones(length(essential.fenids),1);
        end
        if (~isfield( essential, 'component' )) || isempty(essential.component)
            essential.component =1:model_data.un1.dim;
        end
        fixed_value =0;
        if (isfield( essential, 'fixed_value' ))
            fixed_value = essential.fixed_value;
        end
        val=zeros(length(essential.fenids),1)+fixed_value;
        for k=1:length( essential.component)
            model_data.un1 = set_ebc(model_data.un1, essential.fenids, essential.is_fixed, essential.component(k), val);
        end
        model_data.un1 = apply_ebc (model_data.un1);
        % If the essential boundary condition is to be time-dependent, the
        % values must be supplied by a function.   Find out…
        essential.time_dependent =(isa(essential.fixed_value,'function_handle'));
        model_data.boundary_conditions.essential{j} =essential;
    end
    clear essential fenids fixed component fixed_value  val
end

% Number the equations: The displacements may be prescribed
model_data.un1 = numberdofs (model_data.un1);
% and that would also be reflected in the velocity vector..
model_data.v = model_data.un1;
model_data.un = model_data.un1;
clear un1 un

% Associate geometry and compute the mass matrix
M=  sparse(model_data.un1.nfreedofs,model_data.un1.nfreedofs);
for i=1:length(model_data.region)
    region =model_data.region{i};
    if (~isfield(region, 'femm' ))
        error('Must supply the finite element machine');
    end
    % Give the  FEMM a chance  to precompute  geometry-related quantities
    region.femm = associate_geometry(region.femm,model_data.geom);
    M = M + lumped_mass(region.femm, sysmat_assembler_sparse, model_data.geom, model_data.un1); model_data.region{i}=region;
    clear region Q prop mater Rm  femm
end

%     Set the initial condition
if isfield(model_data,'initial_condition')
    if (isa(model_data.initial_condition.u_fixed_value,'function_handle'))
        model_data.un1.values = model_data.initial_condition.u_fixed_value(model_data.geom.values);
    else
        uval=gather_sysvec(model_data.un1)*0+model_data.initial_condition.u_fixed_value;
        model_data.un1 = scatter_sysvec(model_data.un1,uval);
    end
    if (isa(model_data.initial_condition.v_fixed_value,'function_handle'))
        model_data.v.values = model_data.initial_condition.v_fixed_value(model_data.geom.values);
    else
        vval=gather_sysvec(model_data.v)*0+model_data.initial_condition.v_fixed_value;
        model_data.v = scatter_sysvec(model_data.v,uval);
    end
else % no initial conditions were supplied, assume  all displacements and velocities are zero.
    model_data.un1=gather_sysvec(model_data.un1)*0;
    model_data.v=gather_sysvec(model_data.v)*0;
end

% Process the traction boundary condition
if (isfield(model_data, 'boundary_conditions' ))&&(isfield(model_data.boundary_conditions, 'traction' ))
    for j=1:length(model_data.boundary_conditions.traction)
        traction =model_data.boundary_conditions.traction{j};
        traction.femm = femm_deformation_linear (struct ('material',[],...
            'fes',traction.fes,...
            'integration_rule',traction.integration_rule));
        % If the traction boundary condition is to be time-dependent, the
        % values must be supplied by a function.   Find out!
        traction.time_dependent =(isa(traction.traction,'function_handle'));
        model_data.boundary_conditions.traction{j}=traction;
    end
    clear traction fi  femm
end

% Solve

if isempty(dt)
    % Find the stable time step.  Compute the largest eigenvalue (angular
    % frequency of vibration), and determine the time step from it.   
    dt=inf;
    for i=1:length(model_data.region)
        region =model_data.region{i};
        stabldt = estimate_stable_step (region.femm, model_data.geom+model_data.un1);
        dt =dt_reduction * min([dt, stabldt]);
        clear region
    end
end
model_data.dt=dt;

% Let us begin:
t=0; 
% Initial displacement, velocity, and acceleration.
U0 = gather_sysvec(model_data.un1);
V0= gather_sysvec(model_data.v);
A0 =[];% The acceleration will be computed from the initial loads.
F0 = 0*U0;% Zero out the load
% Process the time-independent traction boundary condition
if (isfield(model_data.boundary_conditions, 'traction' ))
    for j=1:length(model_data.boundary_conditions.traction)
        traction =model_data.boundary_conditions.traction{j};
        if (~traction.time_dependent)
            fi= force_intensity(struct('magn',traction.traction));
            F0 = F0 + distrib_loads(traction.femm, sysvec_assembler, model_data.geom, model_data.un1, fi, 2);
        end
    end
    clear traction fi
end
%         % Process time-dependent essential boundary conditions:
if ( isfield(model_data,'boundary_conditions')) && ...
        (isfield(model_data.boundary_conditions, 'essential' ))
    for j=1:length(model_data.boundary_conditions.essential)
        essential =model_data.boundary_conditions.essential{j};
        if (~essential.time_dependent)% this needs to be computed
            % only for truly time-independent EBC. The field needs to be
            % cleared of fixed values for each specification of the
            % essential boundary conditions in order not to apply them
            % twice.
            fixed_value =essential.fixed_value;
            for k=1:length(essential.component)
                model_data.un1 = set_ebc(model_data.un1, essential.fenids, ...
                    essential.is_fixed, essential.component(k), fixed_value);
            end
            model_data.un1 = apply_ebc (model_data.un1);
            for i=1:length(model_data.region)
                F0 = F0 + nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.un1);
            end
            model_data.un1.fixed_values(:)=0;% Zero out the fixed values for the next load
        end
    end
    clear essential
end
% This is the time-independent load vector
F0indep=F0;
%     First output is the initial condition
if ~isempty(observer)% report the progress
    model_data.dt=dt;
    observer (t,model_data);
end
step =0;
while t <tend
    step = step  +1;
    % If desired, estimate the stable time step
    if (mod(step,steps_between_dt_estimations)==0)
        dt=inf;
        for i=1:length(model_data.region)
            stabldt = estimate_stable_step (model_data.region{i}.femm, model_data.geom+model_data.un1);
            dt =dt_reduction * min([dt, stabldt]);
        end
        dt
    end
    
    F0 = 0*F0;% Zero out the load
    F0 = F0 + F0indep;% Add on the time-independent load vector
    %         % Process time-dependent essential boundary conditions:
    if ( isfield(model_data,'boundary_conditions')) && ...
            (isfield(model_data.boundary_conditions, 'essential' ))
        Any_essential_time_dependent =false;;
        for j=1:length(model_data.boundary_conditions.essential)
            essential =model_data.boundary_conditions.essential{j};
            if (essential.time_dependent)% this needs to be computed
                % only for truly time-dependent EBC
                Any_essential_time_dependent =true;;
                if (isa(essential.fixed_value,'function_handle'))
                    fixed_value =essential.fixed_value(t);
                else
                    fixed_value =essential.fixed_value;
                end
                for k=1:length(essential.component)
                    model_data.un1 = set_ebc(model_data.un1, essential.fenids, ...
                        essential.is_fixed, essential.component(k), fixed_value);
                end
                model_data.un1 = apply_ebc (model_data.un1);
            end
        end
        clear essential
        % Add the loads due to fixed displacements
        if (Any_essential_time_dependent)% Only if any time-dependent
            for i=1:length(model_data.region)
                F0 = F0 + nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.un1);
            end
        end
    end
    % Process the traction boundary condition
    if (isfield(model_data, 'boundary_conditions' ))&&(isfield(model_data.boundary_conditions, 'traction' ))
        for j=1:length(model_data.boundary_conditions.traction)
            traction =model_data.boundary_conditions.traction{j};
            if (traction.time_dependent)
                fi= force_intensity(struct('magn',traction.traction(t)));
                F0 = F0 + distrib_loads(traction.femm, sysvec_assembler, model_data.geom, model_data.un1, fi, 2);
            end
        end
        clear traction fi
    end
    % If this is the first step compute the initial acceleration.
    if (isempty(A0))
        A0=M\F0;
    end
    % Add  time-dependent traction loads [to be done]
    U1 = U0 +dt*V0+(dt^2)/2*A0;% displacement update
    model_data.un1 = scatter_sysvec(model_data.un1, U1);
    % Add the restoring forces
    for i=1:length(model_data.region)
        region =model_data.region{i};
        [FR,region.femm]=restoring_force(region.femm,sysvec_assembler, model_data.geom,model_data.un1,model_data.un,dt);
        F1 = F0 + FR;
        model_data.region{i}=region;
        clear region
    end
    % Compute the new acceleration.
    A1=M\(F1);
    % Update the velocity
    V1 = V0 +(dt/2)* (A0+A1);
    % Bring the the displacement and velocity fields up to date
    model_data.v = scatter_sysvec(model_data.v, V1);
    % Switch the temporary vectors for the next step.
    U0 = U1;
    V0 = V1;
    A0 = A1;
    model_data.un =  model_data.un1;
    if (t==tend)
        break;
    end
    if (t+dt>tend)
        dt =tend-t;
    end
    t=t+dt;
    if ~isempty(observer)% report the progress
        model_data.dt=dt;
        observer (t,model_data);
    end
end

end
