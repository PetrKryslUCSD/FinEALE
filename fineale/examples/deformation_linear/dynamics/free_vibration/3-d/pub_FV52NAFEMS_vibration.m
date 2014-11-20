%% Vibration analysis of simply-supported thick  (solid) plate
%

%%
% Link to the  <matlab:edit('pub_FV52NAFEMS_vibration') m-file>.
%

%% Description
%
% Free-vibration problem is solved for a homogeneous square plate,
% simply-supported on the circumference.
% This is the NAFEMS Benchmark, Test No. FV52.


%%
% The plate is discretized with solid elements. The simple support
% condition is approximated by distributed rollers on the boundary.
% Because only the out of plane displacements are prevented, the structure
% has three rigid body modes in the plane of the plate.

%%
% The nonzero benchmark frequencies are (in hertz): 45.897,   109.44,  109.44, 167.89,
%  193.59,  206.64,  206.64. Note: these are not only for the out of plane 
% mode shapes. The in-plane modes (5,6,7) are included.

%% Solution
%
function pub_FV52NAFEMS_vibration
    
    %%
    % Define the material properties.
    pu=physical_units_struct;
    % Parameters:
    E = 200e3*pu.MEGA*pu.PA;
    nu = 0.3;
    rho= 8000*pu.KG/pu.M^3;
    
    %%
    %
    L =10*pu.M;% span of the plate
    t =1.0*pu.M;% thickness of the plate
    
    %%
    % The chosen mesh parameters. This is the  coarse mesh as specified in the benchmark.
    nL= 4;% number of elements span wise
    nt = 1;% number of elements through the thickness
    
    
    %%
    % The mesh is generated
    [fens,fes] = H8_block(L,L,t,nL,nL,nt);;
    
    %%
    % The chosen elements are the serendipity hexahedra.
    [fens,fes] = H8_to_H20(fens,fes);
    
    
    %%
    % We are ready to bundle up the model data so they can be passed to the solver.
    clear model_data
    model_data.fens =fens;% the finite element node set
    
    %%
    % Note that we are specifying the material parameters and the material
    % orientation matrix. The integration rule is going to be used for both the
    % stiffness matrix and the mass matrix.
    clear region
    region.rho =rho;
    region.E=E;
    region.nu=nu;
    region.fes= fes;% set of finite elements for the interior of the domain
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
    model_data.region{1} =region;
    
    %%
    %
    % The support conditions approximate simply-supported edges.  All the
    % sides of the plate are fixed in the transverse direction  (Z
    % displacement).
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
        'inflate',0.001*t)),fenode_select(fens, struct('box', [L,L,-Inf,Inf,-Inf,Inf],...
        'inflate',0.001*t)),fenode_select(fens, struct('box', [0,L,0,0,-Inf,Inf],...
        'inflate',0.001*t)),fenode_select(fens, struct('box', [0,L,L,L,-Inf,Inf],...
        'inflate',0.001*t))];
    model_data.boundary_conditions.essential{1} = essential;
    
    %%
    % How many natural frequencies should be calculated?
    model_data.neigvs= 10;
    %%
    % Three rigid body modes are to be expected with the present boundary
    % conditions.  We have to use mass shifting. 10 Hz appears to be a good
    % frequency between the first nonzero natural frequency  and the rigid body
    % mode frequency of zero.
    model_data.omega_shift= 10*2*pi;
    
    %%
    % The modal analysis solver is now ready to be invoked.
    
    model_data = deformation_linear_modal_analysis(model_data);
    
    
    %%
    % Due to the mass-shifting,  the frequencies may come out  with  nonzero
    % (but hopefully small) imaginary parts.  Remove the imaginary parts as
    % they have no meaning.
    format short e
    model_data.Omega'
    model_data.Omega  = real(model_data.Omega);
    model_data.W  = real(model_data.W);
    
 
%% 
% Furthermore, let us get rid of the rigid body modes (first three).
model_data.Omega  = model_data.Omega(4:end);
    model_data.W  = model_data.W(:,4:end);
    
    %%
    % The modal-plot algorithm can be called to produce the plot of the
    % second and third natural  frequency mode.
    model_data.postprocessing.u_scale= 2;
    model_data.postprocessing.modelist=3;
    model_data.postprocessing.cmap=parula;
    model_data=deformation_plot_modes(model_data);
    
    %% Discussion
    %
    %%
    % The fundamental frequency is extracted from the updated |model_data| data
    % structure. 
    f=model_data.Omega'/2/pi;
    disp(['Natural frequencies '  ': ' num2str(f) ' [Hz]' ]);
    
    %%
    % These computed frequencies should be compared with the benchmark values
    % of 45.897,   109.44,  109.44, 167.89,  193.59,  206.64,  206.64 [Hz].
    % In percent we have the comparison
    format short
    f./[45.897,   109.44,  109.44, 167.89,  193.59,  206.64,  206.64]*100
end