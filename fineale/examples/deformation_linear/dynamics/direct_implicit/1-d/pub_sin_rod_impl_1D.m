%% Suddenly applied harmonic load on a bar; 1-D model with trapezoidal rule
%

%%
% Link to the  <matlab:edit('pub_sin_rod_impl_1D') m-file>.
%

%% Description
% The structure is a clamped bar that consists of two materials, a short section at the free
% end of stiff material,  and a long section of flexible material adjacent
% to the clamped end (see Figure 1).
%%
% The free end is loaded by a suddenly-applied harmonically varying
% concentrated force of  4.0 N. The deflection of the tip is sought at the time 0.01
% seconds.  The deflection of 0.3558 mm is the reference solution described
% by Bathe, "Finite Element Procedures," Prentice Hall, Englewood Cliffs,
% New Jersey, 1996, pp.818-821.


%% 
% The equation of motion is to be integrated with a given time step (0.0004
% seconds) and therefore  we will use an implicit (trapezoidal) time
% integration rule.

%%
%
% <html>
% <table border=0><tr><td>
% <img src="../docs/pub_sin_rod_1D.jpg" width = "60%">
% </td></tr>
% <tr><td>Figure 1. Definition of the geometry of the two-material rod with a harmonic force suddenly applied at the free tip</td></tr>
% </table>
% </html>


%% Solution
%
function pub_sin_rod_impl_1D
    u=physical_units_struct;
    %%
    % Define the material properties of the stiff and flexible sections.  Note
    % that even though we are defining nonzero Poisson's ratios, the solution
    % is going to be obtained with a uniaxial model.  The uniaxial material
    % model extracts from the 3-D material stiffness matrix only the Young's
    % modulus by the model reduction procedure based on the assumption that the
    % transverse stretches are nonzero but the transverse normal stresses (and
    % all shears) are zero.
    E1=2e5*u.MEGA*u.PA;%  stiff (short) section
    nu1=0.3;
    rho1  = 7800*u.KG/u.M^3;
    E2=4.432*u.MEGA*u.PA;% flexible (long) section of the rod
    nu2=0.3;
    rho2  = 1560*u.KG/u.M^3;
    %%
    % Geometrical dimensions of the rod..
    L = 1.0*u.M;% total length of the rod
    L1= 0.05*u.M;%length of the stiff section of the rod
    W = 0.02*u.M;% cross-sectional dimension
    H = W;
    tolerance  =W/1000;% geometrical tolerance
    
    %%
    %  The target time and the time step are defined.
    tend = 0.01*u.SEC;
    dt=0.0004*u.SEC;
    
    %%
    % The mesh is going to consist of 20 elements.
    n=20;
    xs=[linspace(0,L,n+1),L1]';% locations of the nodes
    %%
    % Make sure the interface between the stiff and flexible section  is a
    % location of a node. Note that for the finite elements representing the
    % volume of the rod we are setting the cross-sectional area (W*H) as the "other
    % dimension".
    xs =unique_within_tolerance(xs,tolerance);
    [fens,fes] = L2_blockx(xs,struct('other_dimension',W*H));
    
    
    %%
    % Compose the model data.
    clear model_data
    model_data.fens =fens;
    
    %%
    % The first, stiff, section of the rod consists of finite elements between
    % 0 and the length of the stiff section.
    clear region
    region.property = 'isotropic';
    region.E =E1;
    region.nu=nu1;
    region.rho=rho1;
    rl1=fe_select (fens,fes,struct('box',[0,L1],'inflate',tolerance));
    region.fes= subset(fes,rl1);
    region.integration_rule = gauss_rule (struct('dim', 1, 'order', 2));
    model_data.region{1} =region;
    
    %%
    % The flexible section of the rod consists of the finite elements that remain.
    clear region
    region.property = 'isotropic';
    region.E =E2;
    region.nu=nu2;
    region.rho=rho2;
    region.fes= subset(fes,setdiff(1:count(fes),rl1));
    region.integration_rule = gauss_rule (struct('dim', 1, 'order', 2));
    model_data.region{2} =region;
    
    %%
    % The essential boundary condition at the clamped end.
    clear essential
    essential.component= [];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[L,L],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    %%
    % The attraction boundary condition at the free end.  Given that the total
    % force is to be 4.0 Newton the traction (force density) is computed from the
    % cross-sectional area.
    clear traction
    traction.traction= @(t)4*u.NT/(W*H)*u.PA*sin(150*t);
    
    %%
    % The other dimension for the finite elements on the boundary (they are the
    % point elements, P1) is the cross-sectional area.
    bfes=mesh_boundary(fes,struct('other_dimension',H*W)) ;
    ll=fe_select (fens,bfes,struct('box',[0,0],'inflate',tolerance));
    traction.fes=subset(bfes,ll);
    traction.integration_rule = point_rule ;
    model_data.boundary_conditions.traction{1} = traction;
    
    %%
    % The rod is initially at rest, no deformation.
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)0*x;
    model_data.initial_condition = initial_condition;
    
    %%
    % Finally we set the control parameters.
    model_data.tend = tend;;
    model_data.dt = dt;;
    %%
    % The observer function collects the displacement at the tip of the rod as
    % a function of time.
    tip=fenode_select (fens,struct('box',[0 0],'inflate',tolerance));
    tipu  = [];
    function output(t, model_data)
        tipu=[tipu,gather_values(model_data.u,tip)];
    end
    model_data.observer  =@output;;
    
%% 
% Now we call the transient trapezoidal-rule solver. Only an implicit
% method can deal with their given time step as it is above the critical
% time step for  explicit methods (as determined by the stiff portion of
% the rod).
    model_data =deformation_linear_direct_implicit_TRAP_Rayleigh(model_data);
    
%% 
% The displacement at the end of the monitored time interval is:
    tipu(end)/u.MM
%% 
% The deflection curve at the end of thetime interval is plotted:
    plot(model_data.geom.values,model_data.u.values);
    set(gca,'xlim',[0,L]);
    set(gca,'ylim',[0,5e-4]*1);
    labels(  'X [m]', 'Deflection [m]')
%% Discussion
% 
end