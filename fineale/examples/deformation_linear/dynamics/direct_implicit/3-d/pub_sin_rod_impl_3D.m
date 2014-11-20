%% Suddenly applied harmonic load on a bar; 3-D model with implicit integration
%

%%
% Link to the  <matlab:edit('pub_sin_rod_impl_3D') m-file>.

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
% The equation of motion will be integrated with the trapezoidal
% time integration rule, and  the time step will be taken as 0.0004 seconds.

%%
% The bar is assumed to have a circular cross-section and the finite
% element mesh is generated to orient the cylinder along the X axis.  Refer
% to Figure 1.
%%
%
% <html>
% <table border=0><tr><td>
% <img src="../docs/pub_sin_rod_3D.jpg" width = "60%">
% </td></tr>
% <tr><td>Figure 1. Definition of the geometry of the two-material rod with
% a harmonic force suddenly applied at the free tip</td></tr>
% </table>
% </html>


%% Solution
%
function pub_sin_rod_impl_3D
    u=physical_units_struct;
    %%
    % Define the material properties of the stiff and flexible sections.
    E1=2e5*u.MEGA*u.PA;%  stiff (short) section
    nu1=0.3;
    rho1  = 7800*u.KG/u.M^3;
    E2=4.432*u.MEGA*u.PA;% flexible (long) section of the rod
    nu2=0.3;
    rho2  = 1560*u.KG/u.M^3;
    %%
    % Geometrical dimensions of the rod.
    L = 1.0*u.M;% total length of the rod
    L1= 0.05*u.M;%length of the stiff section of the rod
    A=0.0004*u.M^2; % cross-sectional area
    R = sqrt(A/pi);% radius of the cross-section
    tolerance  =R/100;% geometrical tolerance
    
    %%
    % Applied axial force (with sinusoidal time variation).
    P=   4.0*u.NT;% total applied force
    
    %%
    %  The target time is defined.
    tend = 0.01*u.SEC;
    
    %%
    %  The time step is defined.
    dt = 0.0004*u.SEC;
    
    %%
    % The mesh is going to consist of 20 elements longitudinally and two elements radially.
    nL=20; nR=2;
    [fens,fes] = H8_cylinder_n(R, L, nR, nL);
    %%
    % The default orientation of the cylinder is to align its axis with the Z
    % direction.  Therefore we need to rotate the mesh by 90° to align the axis
    % of the cylinder with the X-axis.
    [fens] = rotate_mesh(fens, [0;90/180*pi;0], [0;0;0]);;
    
    %%
    % If desired, the mesh may be rendered.
    %       drawmesh({fens,fes},'fes','facecolor','red');  labels
    %%
    % Now we prepare the model data.
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
    rl1=fe_select (fens,fes,struct('box',[0,L1,-inf,inf,-inf,inf],'inflate',tolerance));
    region.fes= subset(fes,rl1);
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
    model_data.region{1} =region;
    
    %%
    % The flexible section of the rod consists of the finite elements that remain.
    clear region
    region.property = 'isotropic';
    region.E =E2;
    region.nu=nu2;
    region.rho=rho2;
    region.fes= subset(fes,setdiff(1:count(fes),rl1));
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
    model_data.region{2} =region;
    
    %%
    % The essential boundary condition at the clamped end. We select all the nodes  near the plane $x=L$.
    clear essential
    essential.component= [];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[L,L,-inf,inf,-inf,inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    %%
    % The traction boundary condition at the free end.  Given that the total
    % force is to be 4.0 Newton the traction (force density) is computed from the
    % cross-sectional area.
    clear traction
    traction.traction= @(t)[P/A*sin(150*t);0;0];
    
    %%
    % The finite elements on the boundary are quadrilaterals.
    bfes=mesh_boundary(fes) ;
    ll=fe_select (fens,bfes,struct('box',[0,0,-inf,inf,-inf,inf],'inflate',tolerance));
    traction.fes=subset(bfes,ll);
    traction.integration_rule = gauss_rule (struct('dim', 2, 'order', 2));
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
    model_data.dt = 0.0004*u.SEC;
    %%
    % The observer function collects the displacement at the tip of the rod as
    % a function of time.
    tip=fenode_select (fens,struct('box',[0,0,-inf,inf,-inf,inf],'inflate',tolerance));
    tipu  = []; ts=[];
    function output(t, model_data)
        cutip=gather_values(model_data.u,tip);
        tipu=[tipu,cutip(:,1)]; ts=[ts,t];
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
    mean(tipu(:,end))/u.MM
    
    %%
    % which is to be compared with the reference deflection of 0.3558 mm.
    %%
    % The deflection curve at the end of the time interval is plotted:
    plot(model_data.geom.values(:,1),model_data.u.values(:,1),'+', 'linewidth',2);
    set(gca,'xlim',[-eps,L]);
    set(gca,'ylim',[0,4e-4]*1);
    labels(  'Axial location of node [m]', 'Axial displacement at node [m]')
    %% Discussion
    %
    % In addition to the longitudinal discretization and  the time step,
    % the 3-D model is also affected by the discretization of the
    % cross-section.  The relatively crude  mesh in the cross-section
    % affects the total applied force, and therefore the deflection is
    % underestimated.  This error would be reduced by using more elements
    % through the cross-section.
end