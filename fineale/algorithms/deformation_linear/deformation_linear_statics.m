function model_data = deformation_linear_statics(model_data)
% Algorithm for static linear deformation (stress) analysis.
%
% function model_data = deformation_linear_statics(model_data)
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
% For multi point constraints (MPC) (optional):
% model_data.mpc= cell array of structs, each for one MPC.
%      node_list = list of node numbers involved in the MPC,
%      dof_list= numbers of degrees of freedom for the nodes above,
%      umultipliers=multipliers for the nodes above, 
%      penfact=the penalty factor to multiply  the constraint matrix,
%          The MPC looks like this: sum_i m_i u_{dof(i),node(i)} =0
%          where m_i is the multiplier.
%
% 
% Control parameters:
% The following attributes  may be supplied as fields of the model_data struct:
%      renumber = true or false flag (default is true)
%      renumbering_method = optionally choose the renumbering
%           method  ('symrcm' or 'symamd')
%      factorize = should the solution method be based on the
%           Choleski factorization? true or false flag (default is true)
% 
% Output
% model_data =  structure that was given on input updated with:
% The struct model_data on output incorporates in addition to the input the fields
%     geom = geometry field
%     u = computed displacement field
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
    
    factorize = true; % Should we factorize?
    if (isfield(model_data,'factorize'))
        factorize  =model_data.factorize;
    end
    
    
    % Extract the nodes
    fens =model_data.fens;
    
    % Construct the geometry field
    geom = nodal_field(struct('name',['geom'], 'dim', fens.dim, 'fens', fens));
    
    % Construct the displacement field
    u = nodal_field(struct('name',['u'], 'dim', geom.dim, 'nfens', geom.nfens));
    
    % Apply the essential boundary conditions on the displacement field
    if (isfield(model_data.boundary_conditions, 'essential' ))
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
            fixed_value =0;
            if (isfield( essential, 'fixed_value' ))
                if (strcmp(class(essential.fixed_value),'function_handle'))
                   fixed_value = essential.fixed_value(geom.values(fenids,:));
                else
                    fixed_value = essential.fixed_value;
                end
            end
            if (length(fixed_value)==1), fixed_value =repmat(fixed_value,length(fenids),1); end
            val=zeros(length(fenids),1)+fixed_value;
            for k=1:length( component)
                u = set_ebc(u, fenids, is_fixed, component(k), val);
            end
            u = apply_ebc (u);
        end
        clear essential fenids fixed component fixed_value  val
    end
    
    % Number the equations
    u = numberdofs (u, Renumbering_options);
    
    % Initialize the heat loads vector
    F =zeros(u.nfreedofs,1);
    
    % Construct the system stiffness matrix
    K=  sparse(u.nfreedofs,u.nfreedofs);
    for i=1:length(model_data.region)
        region =model_data.region{i};
        % Create  the finite element model machine, if not supplied
        if (~isfield( region, 'femm'))
            if (isfield(region, 'property' ))
                if (strcmp( region.  property, 'orthotropic' ))
                    prop = property_deformation_linear_ortho (...
                        struct('E1',region.E1,'E2',region.E2,'E3',region.E3,'G12',region.G12,'G13',region.G13,'G23',region.G23,'nu12',region.nu12,'nu13',region.nu13,'nu23',region.nu23));
                elseif (strcmp( region.  property, 'isotropic' ))
                    prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
                else% default
                    prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
                end
            else
                prop=property_deformation_linear_iso(struct('E',region.E,'nu',region.nu));
            end
            if (u.dim==1)
                mater = material_deformation_linear_uniax (struct('property',prop ));
            elseif (u.dim==2)
                mater = material_deformation_linear_biax (struct('property',prop, ...
                    'reduction',region.reduction));
            else
                mater = material_deformation_linear_triax (struct('property',prop ));
            end
            Rm = [];%  Is the material orientation matrix supplied?
            if (isfield( region, 'Rm'))
                Rm= region.Rm;
            end
            region.femm = femm_deformation_linear (struct ('material',mater,...
                'fes',region.fes,...
                'integration_rule',region.integration_rule,'Rm',Rm));
        end
        % Give the  FEMM a chance  to precompute  geometry-related quantities
        region.femm = associate_geometry(region.femm,geom);
        % Compute and assemble the stiffness matrix.
        K = K + stiffness(region.femm, sysmat_assembler_sparse, geom, u); 
        % Loads due to the essential boundary conditions on the displacement field
        if (isfield(model_data.boundary_conditions, 'essential' ))
            F = F + nz_ebc_loads(region.femm, sysvec_assembler, geom, u);
        end
        model_data.region{i}=region;
        clear region Q prop mater Rm  femm
    end
     
    % Process the body load
    if (isfield(model_data, 'body_load' ))
        for j=1:length(model_data.body_load)
            body_load =model_data.body_load{j};
            femm = femm_deformation_linear (struct ('material',[],...
                'fes',body_load.fes,...
                'integration_rule',body_load.integration_rule));
            fi= force_intensity(struct('magn',body_load.force));
            F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 3);        
        end
        clear body_load fi  femm
    end
    
    % Process the traction boundary condition
    if (isfield(model_data.boundary_conditions, 'traction' ))
        for j=1:length(model_data.boundary_conditions.traction)
            traction =model_data.boundary_conditions.traction{j};
            femm = femm_deformation_linear (struct ('material',[],...
                'fes',traction.fes,...
                'integration_rule',traction.integration_rule));
            fi= force_intensity(struct('magn',traction.traction));
            F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 2);        
        end
        clear traction fi  femm
    end
    
    % Process the nodal force boundary condition
    if (isfield(model_data.boundary_conditions, 'nodal_force' ))
        for j=1:length(model_data.boundary_conditions.nodal_force)
            nodal_force =model_data.boundary_conditions.nodal_force{j};
            femm = femm_deformation_linear (struct ('material',[],...
                'fes',fe_set_P1(struct('conn',reshape(nodal_force.node_list,[],1))),...
                'integration_rule',point_rule));
            fi= force_intensity(struct('magn',nodal_force.force));
            F = F + distrib_loads(femm, sysvec_assembler, geom, u, fi, 0);          
        end
        clear nodal_force fi femm
    end
    
    % Apply multi point constraints
    if isfield(model_data,'mpc')
        for i=1:length(model_data.mpc)
            mpc =model_data.mpc{i};
            dofnums=0*mpc.umultipliers;% Construct an array of the degree of freedom numbers
            for kx=1:length(mpc.node_list)
                dofnums(kx)=u.dofnums(mpc.node_list(kx),mpc.dof_list(kx));
            end
            % Now call the utility function to calculate the constraint matrix
            [Kmpc,Fmpc]=apply_penalty_mpc(u.nfreedofs,dofnums,mpc.umultipliers,0.0,mpc.penfact);
            K = K + Kmpc;
            F = F + Fmpc;
        end
        clear Kmpc Fmpc
    end
    
    % Solve the system of linear algebraic equations
    if (factorize)
        % Factorize the left-hand side matrix (Choleski)
        K = (K + K')/2;
        [K,status,prm] = chol(K,'lower','vector');%
        if ( status ~= 0 ) error('Choleski factorization failed'), end
        Have_CXSparse = ~ (strcmp('',which ('cs_lsolve')));
        usol=F;
        if (Have_CXSparse)
            usol(prm) = cs_ltsolve (K, cs_lsolve (K, F(prm))); 
        else
            usol(prm)=(K')\(K\F(prm));
        end
    else
        usol= K\F;
    end
    u = scatter_sysvec(u, usol);
    
    % Update the model data
    model_data.geom = geom;
    model_data.u = u;
    model_data.work =dot(F,usol)/2;
    
end
