function model_data = deformation_linear_modal_analysis(model_data)
% Modal (free-vibration) analysis solver.
%
% function model_data = deformation_linear_modal_analysis(model_data)
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
%     integration_rule =integration rule; alternatively, one could 
%           specify separate integration rules for the stiffness and 
%           for the mass matrix.
%     integration_rule_stiffness =integration rule for the stiffness
%     integration_rule_mass =integration rule for the mass
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
%     fixed_value=fixed (prescribed) displacement (scalar)
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
% For multi point constraints (MPC) (optional):
% model_data.mpc= cell array of structs, each for one MPC.
%      mpc.node_list = list of node numbers involved in the MPC,
%      mpc.dof_list= numbers of degrees of freedom for the nodes above,
%      mpc.umultipliers=multipliers for the nodes above, 
%      mpc.penfact=the penalty factor to multiply  the constraint matrix,
%          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
%          where m_i is the multiplier.
%
% Control parameters:
% model_data.neigvs = number of eigenvalues/eigenvectors to compute
% model_data.omega_shift= angular frequency shift for mass shifting
% model_data.renumber = true or false flag (default is true)
% model_data.renumbering_method = optionally choose the renumbering 
%       method  (symrcm or symamd)
% model_data.use_lumped_mass = true or false?  (Default is false: consistent
%       mass)
%
% Output
% model_data = updated structure that was given on input
% The struct model_data on output incorporates in addition to the input the fields
%     model_data.geom = geometry field
%     model_data.u = displacement field
%     model_data.neigvs=Number of computed eigenvectors 
%     model_data.W = Computed eigenvectors, neigvs columns 
%     model_data.Omega=  Computed angular frequencies, array of length neigvs

    
    %     Control parameters
    neigvs=2;
    if (isfield(model_data,'neigvs'))
        neigvs=model_data.neigvs;
    end
    
    omega_shift=0;
    if (isfield(model_data,'omega_shift'))
        omega_shift=model_data.omega_shift;
    end
    
    renumber = true; % Should we renumber?
    if (isfield(model_data,'renumber'))
        renumber  =model_data.renumber;
    end
    Renumbering_options =struct( [] );
    
    use_factorization =true;% quite effective, especially for very large models
    if (isfield(model_data,'use_factorization'))
        use_factorization  =model_data.use_factorization;
    end
    
    use_lumped_mass =true;% Lumped or consistent mass?
    if (isfield(model_data,'use_lumped_mass'))
        use_lumped_mass  =model_data.use_lumped_mass;
    end
    
    % Should we renumber the nodes to minimize the cost of the solution of
    % the coupled linear algebraic equations?
    if (renumber) && (use_factorization)
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
    
    % Construct the displacement field
    u = nodal_field(struct('name',['u'], 'dim', geom.dim, 'nfens', geom.nfens));
    
    % Apply the essential boundary conditions on the displacement field
    if isfield(model_data,'boundary_conditions')...
            && (isfield(model_data.boundary_conditions, 'essential' ))
        for j=1:length(model_data.boundary_conditions.essential)
            essential =model_data.boundary_conditions.essential{j};
            if (isfield( essential, 'fes' ))
                fenids= connected_nodes(essential.fes);
            else
                fenids= essential.node_list;
            end
            is_fixed=ones(length(fenids),1);
            if (isfield( essential, 'is_fixed' ))
                is_fixed= essential.is_fixed;
            end
            component=[];% Implies all components
            if (isfield( essential, 'component' ))
                component= essential.component;
            end
            if (isempty(component))
                component =1:u.dim;
            end
            val=zeros(length(fenids),1);% All supports are immovable
            for k=1:length( component)
                u = set_ebc(u, fenids, is_fixed, component(k), val);
            end
            u = apply_ebc (u);
        end
        clear essential fenids fixed component fixed_value  val
    end
    
    % Number the equations
    u = numberdofs (u, Renumbering_options);
    
    % Construct the system stiffness and mass matrix
    K=  sparse(u.nfreedofs,u.nfreedofs);
    M=  sparse(u.nfreedofs,u.nfreedofs);
    for i=1:length(model_data.region)
        region =model_data.region{i};
        % Create  the finite element model machine, if not supplied
        if (~isfield( region, 'femm'))
            if (isfield(region, 'property' ))
                if (strcmp( region.  property, 'orthotropic' ))
                    prop = property_deformation_linear_ortho (...
                        struct('E1',region.E1,'E2',region.E2,'E3',region.E3,...
                        'G12',region.G12,'G13',region.G13,'G23',region.G23,...
                        'nu12',region.nu12,'nu13',region.nu13,'nu23',region.nu23));
                elseif (strcmp( region.  property, 'isotropic' ))
                    prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
                else% default
                    prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
                end
            else
                prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
            end
            prop.rho=region.rho;
            if (u.dim==1)
                mater = material_deformation_linear_uniax (struct('property',prop ));
            elseif (u.dim==2)
                mater = material_deformation_linear_biax (struct('property',prop, ...
                    'reduction',region.reduction));
            else
                mater = material_deformation_linear_triax (struct('property',prop ));
            end
            Rm = [];
            if (isfield( region, 'Rm'))
                Rm= region.Rm;
            end
            if (isfield(region, 'integration_rule'))
                region.femm = femm_deformation_linear (struct ('material',mater,...
                    'fes',region.fes,...
                    'integration_rule',region.integration_rule,'Rm',Rm));
            else    
                region.femm_stiffness = femm_deformation_linear (struct ('material',mater,...
                    'fes',region.fes,...
                    'integration_rule',region.integration_rule_stiffness,'Rm',Rm));
                region.femm_mass = femm_deformation_linear (struct ('material',mater,...
                    'fes',region.fes,...
                    'integration_rule',region.integration_rule_mass,'Rm',Rm));    
            end       
        end
        if isfield(region,'femm')
            % Give the  FEMM a chance  to precompute  geometry-related quantities
            region.femm = associate_geometry(region.femm,geom);
            K = K + stiffness(region.femm, sysmat_assembler_sparse, geom, u);
            if (use_lumped_mass)
                M = M + lumped_mass(region.femm, sysmat_assembler_sparse, geom, u);
            else
                M = M + mass(region.femm, sysmat_assembler_sparse, geom, u);
            end
        else
            % Give the  FEMM a chance  to precompute  geometry-related quantities
            region.femm_stiffness = associate_geometry(region.femm,geom);
            K = K + stiffness(region.femm_stiffness, sysmat_assembler_sparse, geom, u);
            if (use_lumped_mass)
                M = M + lumped_mass(region.femm_mass, sysmat_assembler_sparse, geom, u);
            else
                M = M + mass(region.femm_mass, sysmat_assembler_sparse, geom, u);
            end
        end
        model_data.region{i}=region;
        clear region Q prop mater Rm  femm
    end
    
    % Apply multi point constraints
    if isfield(model_data,'mpc')
        for i=1:length(model_data.mpc)
            mpc =model_data.mpc{i};
            dofnums=0*mpc.umultipliers;% Construct an array of the degree of freedom numbers
            for kx=1:length(mpc.node_list)
                dofnums(kx)=u.dofnums(mpc.node_list(kx),mpc.dof_list(kx));
            end
            % Now called the utility function to calculate the constraint
            % matrix
            [Kmpc,Fmpc]=apply_penalty_mpc(u.nfreedofs,dofnums,mpc.umultipliers,mpc.penfact);;
            K =K+Kmpc;
        end
    end
    
    % Options for the eigenproblem solution
    clear evopts
    evopts.p= 2*neigvs;
    evopts.issym= true;
    evopts.tol=eps;
    evopts.maxit= 500;
    evopts.disp= 0;
    
    % Solve
    if (~ use_factorization )
        % This is one way of solving the eigenvalue problem, just pass the matrices
        [W,Omega]= eigs(K+omega_shift*M, M, neigvs, 'SM', evopts);
    else
        % This form uses the factorized matrix and has the potential of being much faster
        % Factorize the left-hand side matrix for efficiency (Choleski)
                [mA,status] = chol(K+omega_shift*M,'lower');%,'vector',prm
                if ( status ~= 0 ) error('Choleski factorization failed'), end
                clear K; % Not needed anymore
                mAt= mA';
                [W,Omega]= eigs(@(bv)mAt\(mA\bv), u.nfreedofs, M, neigvs, 'SM', evopts);
%          [W,Omega]= eig(full(K+omega_shift*M), full(M));
    end
    
    %    Subtract the mass-shifting Angular frequency
    Omega =diag(diag(Omega)-omega_shift);
    if (sum(diag(Omega) < 0))
        warning(['Some (' num2str((sum(diag(Omega) < 0))) ') negative angular frequencies detected']);
        Omega=abs(Omega);
    end
    if (sum(abs(imag(diag(Omega))) ~= 0))
        warning(['Some (' num2str(sum(abs(imag(diag(Omega))) ~= 0)) ') complex angular frequencies detected']);
        Omega=real(Omega);
    end
    %    Sort  the angular frequencies by magnitude.  Make sure all
    %    imaginary parts of the eigenvalues are removed.
    [Ignore,ix] =sort (real(diag(Omega)));
    
    % Update the model data: store geometry
    model_data.geom = geom;
    % Store the displacement field
    model_data.u = u;
    % Number of computed eigenvectors
    model_data.neigvs=neigvs;
    %  Computed eigenvectors: we are ignoring the imaginary part here
    %  because the modal analysis is presumed to have been performed for
    %  an undamped structure
    model_data.W = real(W(:,ix));
    %  Computed angular frequencies
    model_data.Omega=diag(sqrt(Omega(ix,ix)));
    
end
