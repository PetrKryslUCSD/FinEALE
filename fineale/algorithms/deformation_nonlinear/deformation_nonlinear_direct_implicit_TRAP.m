function model_data = deformation_nonlinear_direct_implicit_TRAP(model_data)
% Direct-integration nonlinear implicit solver for mechanical vibration.
%
% Algorithm for dynamic nonlinear deformation with implicit  (trapezoidal rule)
% direct integration.
%
% function model_data = deformation_nonlinear_direct_implicit_TRAP(model_data)
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
%      Rayleigh_stiffness= multiplier of the stiffness matrix to obtain
%           the stiffness-proportional damping matrix. For a given loss factor
%           at a certain  frequency, the stiffness-proportional damping
%           coefficient may be estimated as
%                Rayleigh_stiffness = 2*loss_tangent/(2*pi*frequency);
%      Rayleigh_mass= multiplier of the mass matrix to obtain
%           the mass-proportional damping matrix. For a given loss factor
%           at a certain  frequency, the mass-proportional damping
%           coefficient may be estimated as
%                Rayleigh_mass = 2*loss_tangent*(2*pi*frequency);
%
% Output
% model_data = updated structure that was given on input updated as:
% It incorporates in addition to the input the nodal_fields
%     geom = geometry field
%     u = displacement field


%     Control parameters
u0 =@(t,model_data) disp(['Time ' num2str(t)]);
if ( isfield(model_data,'u0'))
    u0  =model_data.u0;;
end
tend =0;
if ( isfield(model_data,'tend'))
    tend  =model_data.tend;;
end
maxdu_tol=0;% Tolerance on the magnitude  of the largest incremental displacement component
if (isfield(model_data,'maxdu_tol'))
    maxdu_tol  =model_data.maxdu_tol;
end
maxbal_tol=0;% Tolerance on the magnitude  of the out-of-balance force
if (isfield(model_data,'maxbal_tol'))
    maxbal_tol  =model_data.maxbal_tol;
end
Newmarkg=1/2;% Newmark gamma
if (isfield(model_data,'Newmarkg'))
    Newmarkg  =model_data.Newmarkg;
    Newmarkb = 1/4*(1/2+Newmarkg)^2;
end
Newmarkb = 1/4*(1/2+Newmarkg)^2;;% Newmark beta
if (isfield(model_data,'Newmarkb'))
    Newmarkb  =model_data.Newmarkb;
end
iteration_observer =[];%  do nothing
if ( isfield(model_data,'iteration_observer'))
    iteration_observer  =model_data.iteration_observer;;
end
increment_observer =@(lambda,model_data) disp(['Time = ' num2str(lambda)]);
if ( isfield(model_data,'increment_observer'))
    increment_observer  =model_data.increment_observer;;
end

% Extract the nodes
fens =model_data.fens;

% Construct the geometry field
model_data.geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));

% Construct the displacement field
u0 = nodal_field(struct('name',['u'], 'dim', model_data.geom.dim, 'nfens', model_data.geom.nfens));

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
            essential.component =1:u0.dim;
        end
        fixed_value =0;
        % If the essential boundary condition is not  supplied by a
        % function, set it equal to zero initially.
        if (isfield( essential, 'fixed_value' )) && (~isa(essential.fixed_value,'function_handle'))
            fixed_value = essential.fixed_value;
        end
        val=zeros(length(essential.fenids),1)+fixed_value;
        for k=1:length(essential.component)
            u0 = set_ebc(u0, essential.fenids, essential.is_fixed, essential.component(k), val);
        end
        u0 = apply_ebc (u0);
        % If the essential boundary condition is to be time-dependent, the
        % values must be supplied by a function.   Find out.
        essential.time_dependent =(isa(essential.fixed_value,'function_handle'));
        model_data.boundary_conditions.essential{j} =essential;
    end
    clear essential fenids fixed component fixed_value  val
end

% Number the equations: The displacements may be fixed...
u0 = numberdofs(u0);
% and that would also be reflected in the velocity vector...
v0 = 0*u0;
% ... and the acceleration vector.
a0 = 0*u0;
% Incremental displacement field
dchi =clone(u0,'dchi');

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
    region.femm = associate_geometry(region.femm,model_data.geom);
    model_data.region{i}=region;
    clear region prop mater Rm  femm
end

% Construct the system mass matrix
M=  sparse(u0.nfreedofs,u0.nfreedofs);
for i=1:length(model_data.region)
    region =model_data.region{i};
    M = M + mass(region.femm, sysmat_assembler_sparse, model_data.geom, u0);
    model_data.region{i}=region;
    clear region Q prop mater Rm  femm
end


% Process the traction boundary condition
if (isfield(model_data.boundary_conditions, 'traction' ))
    for j=1:length(model_data.boundary_conditions.traction)
        traction =model_data.boundary_conditions.traction{j};
        traction.femm = femm_deformation_linear (struct ('material',[],...
            'fes',traction.fes,...
            'integration_rule',traction.integration_rule));
        % If the traction boundary condition is to be time-dependent, the
        % values must be supplied by a function.   Find out...
        traction.time_dependent =(isa(traction.traction,'function_handle'));
        model_data.boundary_conditions.traction{j}=traction;
    end
    clear traction fi  femm
end


%     Set the initial condition
if isfield(model_data,'initial_condition')
    if (isa(model_data.initial_condition.u_fixed_value,'function_handle'))
        u0.values = model_data.initial_condition.u_fixed_value(model_data.geom.values);
    else
        uval=gather_sysvec(u0)*0+model_data.initial_condition.u_fixed_value;
        u0 = scatter_sysvec(u0,uval);
    end
    if (isa(model_data.initial_condition.v_fixed_value,'function_handle'))
        v0.values = model_data.initial_condition.v_fixed_value(model_data.geom.values);
    else
        vval=gather_sysvec(v0)*0+model_data.initial_condition.v_fixed_value;
        v0 = scatter_sysvec(v0,uval);
    end
else % no initial conditions were supplied, assume  all displacements and velocities are zero.
    u0=0*u0;
    v0=0*v0;
end

% Solve
% Let us begin:
t=0; dt=model_data.dt;
%     First output is the initial condition
model_data.un = u0;
model_data.un1 = u0;
model_data.vn = v0;
model_data.vn1 = v0;
if ~isempty(increment_observer)% report the progress
    increment_observer (t,model_data);
end
step =0;
while t <tend
    t=t+dt;
    step = step  +1;
    % Initialization
    dchi = apply_ebc(dchi);% Apply boundary conditions
    u1 = u0; % guess
    stepdchi = 0*apply_ebc(dchi);% Total increment in current step
    a1 = -(1/Newmarkb/dt)*v0 -(1/2-Newmarkb)/Newmarkb*a0;
    v1 = v0 + dt*((1-Newmarkg)*a0 + Newmarkg*a1);
    %     Auxiliary vectors
    dchipv =gather_sysvec(dt*v0 +(dt^2/2*(1-2*Newmarkb))*a0);
    vpv =gather_sysvec(v0 +(dt*(1-Newmarkg))*a0);
    
    %     %         % Process time-dependent essential boundary conditions:
    %     if ( isfield(model_data,'boundary_conditions')) && ...
    %             (isfield(model_data.boundary_conditions, 'essential' ))
    %         Any_essential_time_dependent =false;;
    %         for j=1:length(model_data.boundary_conditions.essential)
    %             essential =model_data.boundary_conditions.essential{j};
    %             if (essential.time_dependent)% this needs to be computed
    %                 % only for truly time-dependent EBC
    %                 Any_essential_time_dependent =true;;
    %                 if (isa(essential.fixed_value,'function_handle'))
    %                     fixed_value =essential.fixed_value(t);
    %                 else
    %                     fixed_value =essential.fixed_value;
    %                 end
    %                 for k=1:length(essential.component)
    %                     model_data.u = set_ebc(model_data.u, essential.fenids, ...
    %                         essential.is_fixed, essential.component(k), fixed_value);
    %                 end
    %                 model_data.u = apply_ebc (model_data.u);
    %             end
    %         end
    %         clear essential
    %     end
    %     % Process the traction boundary condition
    %     if (isfield(model_data.boundary_conditions, 'traction' ))
    %         for j=1:length(model_data.boundary_conditions.traction)
    %             traction =model_data.boundary_conditions.traction{j};
    %             if (traction.time_dependent)
    %                 fi= force_intensity(struct('magn',traction.traction(t)));
    %                 F1 = F1 + distrib_loads(traction.femm, sysvec_assembler, model_data.geom, model_data.u, fi, 2);
    %             end
    %         end
    %         clear traction fi
    %     end
    
    
    iter=1;
    while 1
        %         This will be the total unbalanced force
        F1 = zeros(dchi.nfreedofs,1);
        
        %         F = loads(nl, sysvec_assembler, dchi);
        
        % Add the restoring forces
        for i=1:length(model_data.region)
            region =model_data.region{i};
            FR=restoring_force(region.femm,sysvec_assembler, model_data.geom,u1,u0,dt);
            F1 = F1 + FR;
            model_data.region{i}=region;
            clear region
        end
        
        
        % Construct the system stiffness matrix
        K=  sparse(u1.nfreedofs,u1.nfreedofs);
        for i=1:length(model_data.region)
            K = K + stiffness(model_data.region{i}.femm, sysmat_assembler_sparse, model_data.geom, u1, u0, dt);
            K = K + stiffness_geo(model_data.region{i}.femm, sysmat_assembler_sparse, model_data.geom, u1, u0, dt);
        end
         
        % The total  unbalanced force  with inertial force contribution
        rhs=F1 +M*((-1/(Newmarkb*dt^2))*gather_sysvec(stepdchi)+(1/(Newmarkb*dt^2))*dchipv);
        
        %         Solve for the displacement increment
        dchi = scatter_sysvec(dchi, (K+(1/(Newmarkb*dt^2))*M)\rhs); % Displacement increment
        
        u1 = u1 + dchi;   % increment displacement
        stepdchi = stepdchi + dchi;
        v1 = v1 + (Newmarkg/Newmarkb/dt)*dchi;
        a1 = a1 + (1/Newmarkb/dt^2)*dchi;
        
       
        % Compute the increment data
        ndu=norm(dchi);
        maxdu  =max(max(abs(dchi.values)));
        maxbal=max(abs(rhs));
        
        % Report  iteration results?
        if (~isempty(iteration_observer))
            model_data.un = u0;
            model_data.un1 = u1;
            model_data.vn = v0;
            model_data.vn1 = v1;
            model_data.dt = dt;
            iteration_observer (t,iter,maxdu,maxbal,dchi,model_data);
        end
        
        % Converged?
        if (maxdu <= maxdu_tol) || (maxbal <= maxbal_tol)
            break; % Breakout of the iteration loop
        end;    
        
        iter=iter+1;
    end %  Iteration loop
    
    
    % Converged.
    % Update the FEMMs
    for i=1:length(model_data.region)
        [~,model_data.region{i}.femm]  =...
            restoring_force(model_data.region{i}.femm, sysvec_assembler, ...
            model_data.geom, u1, u0, dt);
    end
    
    % Bring the the displacement and velocity fields up to date
    model_data.un = u0;
    model_data.un1 = u1;
    model_data.vn = v0;
    model_data.vn1 = v1;
            
    % Report results
    if ~isempty(increment_observer)% report the progress
        increment_observer (t,model_data);
    end
    
    % Switch the fields for the next step.
    u0 = u1;       % update the displacement
    v0 = v1;       % update the velocities
    a0 = a1;       % update the accelerations
    
    %     Have we reached the final time?
    if (t==tend)
        break;
    end
    if (t+dt>tend) % Adjust the last time step so that we exactly reach tend
        dt =tend-t;
    end
    
end


end
