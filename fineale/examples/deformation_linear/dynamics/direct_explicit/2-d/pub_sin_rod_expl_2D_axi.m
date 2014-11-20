%% Suddenly applied harmonic load on a bar; 2-D axially symmetric model with centered differences
%

%%
% Link to the  <matlab:edit('pub_sin_rod_expl_2D_axi') m-file>.

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
% The equation of motion will be integrated with the centered-difference
% time integration rule, and  the time step will be taken as 0.0004 seconds.

%%
% The bar is assumed to have a circular cross-section and therefore axial
% symmetry of the entire problem is assumed.  The section of the bar is
% discretized with quadrilaterals.  The assumption in the toolkit is that
% the axis of symmetry is the Y axis, and the radial coordinate is along
% the X axis.  Refer to Figure 1. Therefore the deformation of the bar in
% this solution occurs in the Y direction.
%%
%
% <html>
% <table border=0><tr><td>
% <img src="../docs/pub_sin_rod_2D_axi.jpg" width = "60%">
% </td></tr>
% <tr><td>Figure 1. Definition of the geometry of the two-material rod with
% a harmonic force suddenly applied at the free tip</td></tr>
% </table>
% </html>


%% Solution
%
function pub_sin_rod_expl_2D_axi
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
    % The mesh is going to consist of 20 elements longitudinally and two elements radially.
    nL=20; nR=2;
    ys=[linspace(0,L,nL+1),L1]';% locations of the nodes
    %%
    % Make sure the interface between the stiff and flexible section  is a
    % location of a node. Note that for the finite elements representing
    % the volume of the rod we are setting the flag indicating that the
    % model-reduction procedure is axially symmetric.
    ys =unique_within_tolerance(ys,tolerance);
    xs=linspace(0,R,nR+1);
    [fens,fes] = Q4_blockx(xs,ys,struct('axisymm',true));
    
    
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
    %%
    % Note that we must inform the solver of the model-dimension reduction to
    % be employed for the two-coordinate (axially symmetric) representation.
    region.reduction ='axisymm';
    region.E =E1;
    region.nu=nu1;
    region.rho=rho1;
    rl1=fe_select (fens,fes,struct('box',[-inf,inf,0,L1],'inflate',tolerance));
    region.fes= subset(fes,rl1);
    %%
    % The integration rule is two-dimensional as we integrate in the section of the rod.
    region.integration_rule = gauss_rule (struct('dim', 2, 'order', 2));
    model_data.region{1} =region;
    
    %%
    % The flexible section of the rod consists of the finite elements that
    % are left over from region 1.
    clear region
    region.property = 'isotropic';
    region.reduction ='axisymm';
    region.E =E2;
    region.nu=nu2;
    region.rho=rho2;
    region.fes= subset(fes,setdiff(1:count(fes),rl1));
    region.integration_rule = gauss_rule (struct('dim', 2, 'order', 2));
    model_data.region{2} =region;
    
    %%
    % The essential boundary condition at the clamped end. We select all
    % the nodes  near the plane $x=L$.
    clear essential
    essential.component= [];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[-inf,inf,L,L],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    %%
    % The traction boundary condition at the free end.  Given that the total
    % force is to be 4.0 Newton the traction (force density) is computed from the
    % cross-sectional area.
    clear traction
    traction.traction= @(t)[0;P/A*sin(150*t);];
    
    %%
    % The finite elements on the boundary are line elements, axially
    % symmetric. Note that we are making sure the other dimension of the
    % elements is computed correctly for the surface integrals.
    bfes=mesh_boundary(fes,struct('axisymm', true)) ;
    ll=fe_select (fens,bfes,struct('box',[-inf,inf,0,0],'inflate',tolerance));
    traction.fes=subset(bfes,ll);
    traction.integration_rule = gauss_rule (struct('dim', 1, 'order', 2));
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
    %%
    % The observer function collects the displacement at the tip of the rod as
    % a function of time.
    tip=fenode_select (fens,struct('box',[-inf,inf,0,0],'inflate',tolerance));
    tipu  = []; ts=[];
    function output(t, model_data)
        cutip=gather_values(model_data.u,tip);
        tipu=[tipu,cutip(:,2)]; ts=[ts,t];
    end
    model_data.observer  =@output;;
    
    %%
    % Now we call the transient centered-difference solver. The time step
    % used in the implicit solution of this problem is above the critical
    % time step for  explicit methods (as determined by the stiff portion
    % of the rod). Therefore we let  the solver figure out the stable time step.
    model_data =deformation_linear_direct_explicit_CD(model_data);
    
    
%% 
% The explicit solver  calculated the stable time step and ran therefore with the step of
 model_data.dt
 
%% 
% which means that the solution was obtained in
length(ts)

%% 
%  steps.
    %%
    % The displacement at the end of the monitored time interval is:
    mean(tipu(:,end))/u.MM
    
    %%
    % which is to be compared with the reference deflection of 0.3558 mm.
    %%
    % The deflection curve at the end of the time interval is plotted:
    plot(model_data.geom.values(:,2),model_data.u.values(:,2),'r+', 'linewidth',2);
    set(gca,'xlim',[-eps,L]);
    set(gca,'ylim',[0,4e-4]*1);
    labels(  'Axial location of node [m]', 'Axial displacement at node [m]')
    %% Discussion
    %
    % The axially symmetric 2-D model delivers accuracy comparable to that
    % of  either the 1D or the 3-D model.  The stable time step is quite
    % short because of the stiff part of the rod.  The explicit
    % centered-difference algorithm is therefore quite expensive compared
    % to the implicit trapezoidal rule.
end