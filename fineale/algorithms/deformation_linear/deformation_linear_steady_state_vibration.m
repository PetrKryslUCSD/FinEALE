function model_data = deformation_linear_steady_state_vibration(model_data)
% Steady-state harmonic vibration solver.
%
% Algorithm for dynamic harmonic steady-state linear deformation with
% complex arithmetic. Rayleigh damping may be specified. Damping may also
% be included for fluid-structure interfaces or
% absorbing-boundary conditions.
%
% function model_data = deformation_linear_steady_state_vibration(model_data)
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
%     fixed_value=fixed (prescribed) displacement (scalar); If specified as
%           nonzero, it is assumed that the  nonzero displacement  has a
%           harmonic  time variation, but the distribution of the
%           displacements and the magnitudes do not depend on the
%           frequency. The fixed_value can be also supplied as a function
%           handle: Handle to a function that takes frequency as argument,
%           and returns an appropriate number of fixed displacement values.
%     fes = finite element set on the boundary to which the condition applies
%               or alternatively
%     node_list = list of nodes on the boundary to which the condition applies
%           Only one of fes and node_list needs to be given.
%    
% Example:
%     clear essential
%     essential.component= [1,3];
%     essential.fixed_value= 0.1;
%     essential.node_list = [[fenode_select(fens, struct('box', [0,a,0,0,-Inf,Inf],...
%         'inflate',tolerance))],[fenode_select(fens, struct('box', [0,a,b,b,-Inf,Inf],...
%         'inflate',tolerance))]];;
%     model_data.boundary_conditions.essential{2} = essential;
%
% Example:
%     clear essential
%     essential.component= [1];
%     essential.fixed_value= @(frequency)cos(2*pi*frequency);
%     essential.fes = mesh_boundary(fes);
%     model_data.boundary_conditions.essential{1} = essential;
%
% For traction boundary conditions (optional):
% model_data.boundary_conditions.traction = cell array of struct,
%           each piece of surface with traction boundary condition gets one
%           element of the array with a struct with the attributes
%     traction=traction (vector), or a handle to a function returning a
%           vector.  The function takes frequency as argument, and must
%           return vector of current values of traction, or a handle to a
%           function that returns the vector of current values of traction.
%     fes = finite element set on the boundary to which the condition applies
%     integration_rule= integration rule
%
% Example:
%     clear traction
%     traction.fes =subset(bdry_fes,bclx);
%     traction.traction= @(x) ([sigmaxx(x);sigmaxy(x);0]);
%     traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));
%     model_data.boundary_conditions.traction{1} = traction;
%
% Example:
%     clear traction
%     traction.fes =subset(bdry_fes,bclx);
%     traction_at_x =@(x) ([sigmaxx(x);sigmaxy(x);0]);
%     traction.traction= @(frequency) ( @(x)sin(2*pi*frequency)*traction_at_x(x) );
%     traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));
%     model_data.boundary_conditions.traction{1} = traction;
%    
% For absorbing boundary conditions (optional):
% [[Description needs to be written]]
%    
% Control parameters:
% The following attributes  may be supplied as fields of the model_data struct:
%      frequencies= frequencies at which response should be calculated
%      observer = handle of an observer function (optional)
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
%      renumber = true or false flag (default is true)
%      renumbering_method = optionally choose the renumbering
%           method  ('symrcm' or 'symamd')
%
% Output
% model_data = structure that was given on input updated as follows:
% It incorporates in addition to the input the nodal_fields as attributes
% of the model_data struct:
%     geom = geometry field
%     u = displacement field
%     frequencies  = list of forcing frequencies
%     u_values = cell array of the values of the displacement nodal field,
%          one for each frequency in the list above

%     Control parameters
Rayleigh_stiffness =0;
if ( isfield(model_data,'Rayleigh_stiffness'))
    Rayleigh_stiffness  =model_data.Rayleigh_stiffness;;
end
Rayleigh_mass =0;
if ( isfield(model_data,'Rayleigh_mass'))
    Rayleigh_mass  =model_data.Rayleigh_mass;;
end
observer =@(t,model_data) disp(['Time ' num2str(t)]);
if ( isfield(model_data,'observer'))
    observer  =model_data.observer;;
end
frequencies =[];
if ( isfield(model_data,'frequencies'))
    frequencies  =model_data.frequencies;
end
tend =0;
if ( isfield(model_data,'tend'))
    tend  =model_data.tend;;
end


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
    
% Extract the nodes
fens =model_data.fens;

% Construct the geometry field
model_data.geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));

% Construct the displacement field
model_data.u = nodal_field(struct('name',['u'], 'dim', model_data.geom.dim, 'nfens', model_data.geom.nfens));

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
            essential.component =1:model_data.u.dim;
        end
        fixed_value =0;
        if (isfield( essential, 'fixed_value' ))
            fixed_value = essential.fixed_value;
        end
        % If the essential boundary condition is to be frequency-dependent, the
        % values must be supplied by a function.   Find out...
        essential.frequency_dependent =(isa(fixed_value,'function_handle'));
        val=zeros(length(essential.fenids),1);
        for k=1:length(essential.component)
            model_data.u = set_ebc(model_data.u, essential.fenids, essential.is_fixed, essential.component(k), val);
        end
        model_data.u = apply_ebc (model_data.u);
        model_data.boundary_conditions.essential{j} =essential;
    end
    clear essential fenids fixed component fixed_value  val
end

% Number the equations: The displacements may be fixed... Also apply
% renumbering, if requested above.
model_data.u = numberdofs (model_data.u, Renumbering_options);

% Construct the system stiffness and mass matrix
K=  sparse(model_data.u.nfreedofs,model_data.u.nfreedofs);
M=  sparse(model_data.u.nfreedofs,model_data.u.nfreedofs);
for i=1:length(model_data.region)
    region =model_data.region{i};
    % Create  the finite element model machine, if not supplied
    if (~isfield( region, 'femm'))
        if (isfield(region, 'property' ))
            if (strcmp( region.  property, 'orthotropic' ))
                prop = property_deformation_linear_ortho (...
                    struct('E1',region.E1,'E2',region.E2,'E3',region.E3,...
                    'G12',region.G12,'G13',region.G13,'G23',region.G23,...
                    'nu12',region.nu12,'nu13',region.nu13,'nu23',region.nu23,'rho',region.rho));
            elseif (strcmp( region.  property, 'isotropic' ))
                prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu,'rho',region.rho));
            else% default
                prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu,'rho',region.rho));
            end
        else
            prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
        end
        prop.rho=region.rho;
        if (model_data.u.dim==1)
            mater = material_deformation_linear_uniax (struct('property',prop ));
        elseif (model_data.u.dim==2)
            mater = material_deformation_linear_biax (struct('property',prop, ...
                'reduction',region.reduction));
        else
            mater = material_deformation_linear_triax (struct('property',prop ));
        end
        Rm = [];
        if (isfield( region, 'Rm'))
            Rm= region.Rm;
        end
        region.femm = femm_deformation_linear (struct ('material',mater,...
            'fes',region.fes,...
            'integration_rule',region.integration_rule,'Rm',Rm));
    end
    % Give the FEMM a chance to save computations: the geometry field is
    % not going to change.
    region.femm = region.femm.associate_geometry(model_data.geom);
    % Now we start the computation of the system matrices.
    K = K + stiffness(region.femm, sysmat_assembler_sparse, model_data.geom, model_data.u);
    M = M + mass(region.femm, sysmat_assembler_sparse, model_data.geom, model_data.u);
    model_data.region{i}=region;
    clear region Q prop mater Rm  femm
end
% At this point Rayleigh damping matrix might have been defined as
% C=  Rayleigh_stiffness*K + Rayleigh_mass*M;%  Compute the damping matrix
% For computational efficiency we use a slightly different formula  below.
% The damping matrix  assembled as C includes only the damping ABC.
% C =  sparse(model_data.u.nfreedofs,model_data.u.nfreedofs);

% This is the frequency-independent load vector
F0 = zeros(size(K,1),1);% Zero out the load
F0d = zeros(size(K,1),1);% Zero out the load
F0a = zeros(size(K,1),1);% Zero out the load

% Process the traction boundary condition
if (isfield(model_data.boundary_conditions, 'traction' ))
    for j=1:length(model_data.boundary_conditions.traction)
        traction =model_data.boundary_conditions.traction{j};
        traction.femm = femm_deformation_linear (struct ('material',[],...
            'fes',traction.fes,...
            'integration_rule',traction.integration_rule));
        % If the traction boundary condition is to be frequency-dependent, the
        % values must be supplied by a function.   Find out...
        traction.frequency_dependent =(isa(traction.traction,'function_handle'));
        if (~traction.frequency_dependent)
            fi= force_intensity(struct('magn',traction.traction));
            F0 = F0 + distrib_loads(traction.femm, sysvec_assembler, ...
                model_data.geom, model_data.u, fi, 2);
        end
        model_data.boundary_conditions.traction{j}=traction;
    end
    clear traction fi  femm
end

%  Process frequency-independent essential boundary conditions: For
%  frequency-independent  distribution and magnitudes of prescribed
%  displacements.
if ( isfield(model_data,'boundary_conditions')) && ...
        (isfield(model_data.boundary_conditions, 'essential' ))
    for j=1:length(model_data.boundary_conditions.essential)
        essential =model_data.boundary_conditions.essential{j};
        if (~essential.frequency_dependent)% 
            fixed_value =essential.fixed_value;
            for k=1:length(essential.component)
                model_data.u = set_ebc(model_data.u, essential.fenids, ...
                    essential.is_fixed, essential.component(k), fixed_value);
            end
            model_data.u = apply_ebc (model_data.u);
        end
    end
    clear essential
    % Add the loads due to fixed displacements
    for i=1:length(model_data.region)
        F0d = F0d + nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.u);
        F0a = F0a + nz_ebc_loads_acceleration(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.u);
    end
end
    
% Solve
u_values=cell(1,length(frequencies));
% Loop over all excitation frequencies
for k=1:length(frequencies)
    frequency=frequencies(k);
    omega =2*pi*frequency;
    % Start with the frequency-independent component...
    F1 = F0 + (-omega^2*F0a + 1i*omega*(Rayleigh_stiffness*F0d + Rayleigh_mass*F0a) + F0d);
    %     ...And then frequency-dependent loads
    %  Process frequency-dependent essential boundary conditions:
    if ( isfield(model_data,'boundary_conditions')) && ...
            (isfield(model_data.boundary_conditions, 'essential' ))
        Any_essential_frequency_dependent =false;;
        for j=1:length(model_data.boundary_conditions.essential)
            essential =model_data.boundary_conditions.essential{j};
            if (essential.frequency_dependent)% this needs to be computed
                % only for truly frequency-dependent EBC
                Any_essential_frequency_dependent =true;;
                fixed_value =essential.fixed_value(frequency);
                for ck=1:length(essential.component)
                    model_data.u = set_ebc(model_data.u, essential.fenids, ...
                        essential.is_fixed, essential.component(ck), fixed_value);
                end
                model_data.u = apply_ebc (model_data.u);
            end
        end
        clear essential
        % Add the loads due to fixed displacements
        if (Any_essential_frequency_dependent)% Only if any frequency-dependent
            for i=1:length(model_data.region)
                F1d = nz_ebc_loads(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.u);
                F1a = nz_ebc_loads_acceleration(model_data.region{i}.femm, sysvec_assembler, model_data.geom, model_data.u);
                F1 = F1 + (-omega^2*F1a + 1i*omega*(Rayleigh_stiffness*F1d + Rayleigh_mass*F1a) + F1d);
            end
        end
    end
    % Process the frequency-dependent traction boundary condition
    if (isfield(model_data.boundary_conditions, 'traction' ))
        for j=1:length(model_data.boundary_conditions.traction)
            traction =model_data.boundary_conditions.traction{j};
            if (traction.frequency_dependent)
                fi= force_intensity(struct('magn',traction.traction(frequency)));
                F1 = F1 + distrib_loads(traction.femm, sysvec_assembler, model_data.geom, model_data.u, fi, 2);
            end
        end
        clear traction fi
    end
    % The damping matrix  assembled as C includes only the damping ABC.
    C =  sparse(model_data.u.nfreedofs,model_data.u.nfreedofs);
    % Process the damping ABC (Absorbing Boundary Condition)
    if (isfield(model_data.boundary_conditions, 'damping_abc' ))
        for j=1:length(model_data.boundary_conditions.damping_abc)
            damping_abc =model_data.boundary_conditions.damping_abc{j};
            damping_abc.frequency_dependent =(isa(damping_abc.damping_abc_impedance,'function_handle'));
            damping_abc_impedance =damping_abc.damping_abc_impedance;
            if (damping_abc.frequency_dependent)
                damping_abc_impedance =damping_abc.damping_abc_impedance(frequency);
            end
            damping_abc.femm = femm_deformation_surface_damping (struct ('material',[],...
                'fes',damping_abc.fes,...
                'damping_abc_impedance',damping_abc_impedance,...
                'integration_rule',damping_abc.integration_rule));
            C = C + surface_damper_abc(damping_abc.femm, sysmat_assembler_sparse, model_data.geom, model_data.u);
            model_data.boundary_conditions.damping_abc{j} =damping_abc;
        end
        clear damping_abc
    end
    
    % Solve  for the complex displacement amplitude.
    % Note:  the final expression is equivalent to
    %     U1 = (-omega^2*M + 1i*omega*C + K)\F1;
    % where
    %     C=  Rayleigh_stiffness*K + Rayleigh_mass*M;%  Compute the damping matrix
    U1 = ((-omega^2+1i*omega*Rayleigh_mass)*M + 1i*omega*C + (1+1i*omega*Rayleigh_stiffness)*K)\F1;
    %     [userview, systemview] = memory;      % << add this line
    %     U1 = mldivide( ((-omega^2+1i*omega*Rayleigh_mass)*M + 1i*omega*C + (1+1i*omega*Rayleigh_stiffness)*K), F1 );
    %     if (sum(isnan(U1))>0)%  This shouldn't happen, unfortunately it sometimes does as there seems to be a bug  in the Matlab engine
    %         save('nanDataSet.mat', 'omega', 'Rayleigh_mass', 'M', 'C', 'Rayleigh_stiffness', ...
    %             'K', 'F1', 'U1', 'userview', 'systemview');     % <<  add this line
    %         error('Failure of the complex-linear-equations solver'  )
    %     end
    
    %     The system is unraveled into purely real solution, and the components of the complex solution are then recombined.
%     rU1 = [-omega^2*M+K, -omega*C-omega*Rayleigh_mass*M-omega*Rayleigh_stiffness*K;
%         -omega*C-omega*Rayleigh_mass*M-omega*Rayleigh_stiffness*K, +omega^2*M-K]\[real(F1); -imag(F1)];
%     U1=rU1(1:length(F1))+1i*rU1(length(F1)+1:end);
    if (sum(isnan(U1))>0)%  This shouldn't happen, unfortunately it sometimes does as there seems to be a bug  in the Matlab engine
        save('nanDataSet.mat');     % <<  add this line
        error('Failure of the complex-linear-equations solver'  )
    end
    
    % Update the displacement field
    model_data.u = scatter_sysvec(model_data.u, U1);
    if ~isempty(observer)% report the progress
        observer (frequency,model_data);
    end
    u_values{k}=model_data.u.values; % save the solution
end% Loop over all excitation frequencies

% Provide comprehensive output for visualization
model_data.frequencies=frequencies;
model_data.u_values=u_values;

end
