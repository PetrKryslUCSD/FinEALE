%% Free harmonic vibration of a deep beam.
%

%%
% Link to the  <matlab:edit('pub_test5_free') m-file>.
%

%% Description
%
% This is a test recommended by the National Agency for Finite Element
% Methods and Standards (U.K.): Test 5 from NAFEMS “Selected Benchmarks
% for Forced Vibration,” R0016, March 1993.
%

%%
% The beam is going to be modeled as a three-dimensional solid, with solid
% elements. The beam geometry is given as

%%
% <<../docs/pub_test5_harmonic_1.jpg>>

%%
% The beam is assumed to be simply supported at the ends transversely
% ($u_y=u_z=0$). The cross-section at $x=0$ is also pinned in the axial
% direction and rotation about the axis of the beam is prevented. In the
% cross-section at $x=L$ both the axial rotation and the axial displacement
% are allowed.
%%
% For the solid 3-D beam this seems impossible to achieve with just the
% mechanism of essential boundary conditions.   In this tutorial we will
% apply multi point constraints (MPCs) to construct the pinned conditions.

%%
% The transverse and axial-rotation constrains at $x=0$ is implemented by
% attaching rollers in the direction of the axes perpendicular to the beam
% in both directions on the entire cross-sectional area.
%%
% The axial constraint at $x=0$ is implemented with a multi-point
% constraint of the form

%%
% $\sum_i u_{ix}=0$

%%
% where $i$ are the nodes on the cross-section at $x=0$.
%%
% The transverse constraint at the pin at $x=L$ is implemented with two multi-point
% constraints of the form

%%
% $\sum_i u_{iy}=0$
%%
% and

%%
% $\sum_i u_{iz}=0$

%%
% where $i$ are the nodes on the cross-section at $x=0$.

%% Solution
%

%%
% The code is placed inside a Matlab function in order to be able to define
% the simulation in a nested function.
function pub_test5_free
    
    %%
    % Define the material properties.
    pu=physical_units_struct;
    % Parameters:
    E = 200e3*pu.MEGA*pu.PA;
    nu = 0.3;
    rho= 8000*pu.KG/pu.M^3;
    
    %%
    % Define the geometry and the geometrical tolerance.
    a=2.00*pu.M; b=2.00*pu.M; L= 10*pu.M;
    tolerance =a/1000;
    
    %%
    % The reference solutions for the natural frequencies are:
    f_natural =[ 42.6580   42.6580   71.2610  125.0000  ...
        148.7200  148.7200  213.8900  287.8400   287.8400];
    %%
    % in hertz, and the classification of the modes is
    f_natural_kind={'flexural','flexural','torsional','extensional',...
        'flexural','flexural','torsional','flexural', 'flexural'};
    
    %%
    % Mesh parameters are defined here: number of element edges per dimension.
    na= 2; nb=  2; nL =10;
    
    
    %%
    % These are the material properties for an isotropic homogeneous material.
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    %%
    % The simulation will be carried out for selected element types using this
    % function.  It sets up the model, runs the free vibration analysis, and
    % reports the results.
    function  simulate(description, mf, femmf,  mode_to_show)
        % Create the mesh and initialize the geometry
        [fens,fes]= mf(L,a,b,nL,na,nb);
        %      Shift the geometry so that the axis of the beam is along the
        %      X-axis of the global coordinate system.
        fens = transform_apply(fens,@(x,d)(x-[0,a/2,b/2]), []);
        % Select the boundary faces for the application of the transverse load.
        bfes= mesh_boundary(fes,[]);
        topl =fe_select (fens,bfes,struct('box', [0,L,-Inf,Inf,b/2,b/2],...
            'inflate',tolerance));
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        % Set the finite element model machine
        clear region
        region.femm= femmf(fes);
        model_data.region{1} =region;
        
        % Both cross-sections of the beam are supported by rollers in the
        % transverse direction. In the 3-D solid this may be accomplished
        % by restraining translation in the direction of the support at all
        % nodes in the cross-section.  However, the support at $x=L$ must
        % allow for free axial rotation, and hence there the  restraint
        % must be applied differently.
        clear essential
        essential.component= [2,3];%restrain both cross-sections with rollers
        essential.fixed_value= 0;
        essential.node_list = [...
            fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
            'inflate',tolerance))];
        model_data.boundary_conditions.essential{1} = essential;
        
        % At $x=0$  we need to restrain the axial translation. This will be
        % accomplished with the MPC as
        clear mpc
        mpc.node_list = [fenode_select(fens, ...
            struct('box', [0,0,-Inf,Inf,-Inf,Inf],'inflate',tolerance))];
        mpc.dof_list=1+zeros(size(mpc.node_list));
        mpc.umultipliers=ones(size(mpc.node_list));
        mpc.penfact=1e16;
        model_data.mpc{1} = mpc;
        % At $x=L$ the transverse pin support will be also implemented using
        % MPC.
        clear mpc
        mpc.node_list = [fenode_select(fens, ...
            struct('box', [L,L,-Inf,Inf,-Inf,Inf],'inflate',tolerance))];
        mpc.dof_list=2+zeros(size(mpc.node_list));
        mpc.umultipliers=ones(size(mpc.node_list));
        mpc.penfact=1e16;
        model_data.mpc{2} = mpc;
        clear mpc
        mpc.node_list = [fenode_select(fens, ...
            struct('box', [L,L,-Inf,Inf,-Inf,Inf],'inflate',tolerance))];
        mpc.dof_list=3+zeros(size(mpc.node_list));
        mpc.umultipliers=ones(size(mpc.node_list));
        mpc.penfact=1e16;
        model_data.mpc{3} = mpc;
        
        % Solve
        model_data.neigvs= 20;
        model_data.omega_shift=0;
        model_data.use_factorization= true;
        model_data = deformation_linear_modal_analysis(model_data);
        % These are the natural frequencies in hertz
        f=model_data.Omega(1:9)'/2/pi;
        disp ([description  ': '])
        for j=1:length(f)
            disp(['Natural frequency ' num2str(j) ' = ' ...
                num2str(f(j)) ' [Hz]' ]);
            disp([' f/f_analytical=' num2str(f(j)./f_natural(j)*100)...
                '% (' f_natural_kind{j} ')']);
        end
        % Plot the resulting mode shape(s) requested.
        model_data.postprocessing.u_scale= 2;
        for j1=1:length(mode_to_show)
            j=mode_to_show(j1);
            disp('====================================================')
            disp(['Natural frequency ' num2str(j) ' = ' ...
                num2str(f(j)) ' [Hz]' ]);
            disp([' f/f_analytical=' num2str(f(j)./f_natural(j)*100)...
                '% (' f_natural_kind{j} ')']);
            model_data.postprocessing.modelist= j;
            model_data=deformation_plot_modes(model_data);
            snapnow;
        end
    end
    
    
    %%
    % The 20-node serendipity hexahedron which is  a very accurate element in
    % general yields good results here.
    description ='H20R';
    mf =@H20_block;
    %%
    % Note that we are specifying integration rule for both the stiffness and
    % the mass matrix.  Hence the mass matrix will be singular.
    femmf =@(fes)femm_deformation_linear(struct('fes',fes, ...
        'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    simulate(description, mf, femmf, 1:8);
    
    %%
    % For this material (compressible)  and this type of structure  (massive)
    % the fully-integrated hexahedron would also work well.  We can therefore
    % use the same element with Gauss rule one order higher. The consistent
    % mass matrix will no longer be singular. 
    femmf =@(fes)femm_deformation_linear(struct('fes',fes, ...
        'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    simulate(description, mf, femmf, 3);
    
    %%
    % It is also possible to  use different quadrature rules for stiffness and
    % mass calculations.  The tutorial pub_test5_free_sep explains how to
    % do this.
    
    %%
    % The quadratic tetrahedron, the workhorse of many linear stress analysis
    % packages, produces good-quality results.
    description ='T10';% tetrahedron
    mf =@T10_block;
    femmf =@(fes)femm_deformation_linear(struct('fes',fes,...
        'material',mater,...
        'integration_rule',tet_rule(struct('npts',4))));
    simulate(description, mf, femmf, 3);
    
    %% Discussion
    %
    
    %%
    % To compare vibration results produced by different finite element models
    % needs to be done  both numerically and visually.  Even  though the
    % frequencies may seem close, the mode shapes may be incorrect.
end