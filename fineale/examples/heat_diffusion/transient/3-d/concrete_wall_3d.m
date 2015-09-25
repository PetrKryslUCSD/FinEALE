% Temperature distribution through a concrete wall  based on sudden increase
% of temperature at one side.
function concrete_wall_3d
    u=physical_units_struct;
    kappa=1.81*u.W/u.K/u.M; % W/K/m
    rho = 2350*u.KG/u.M^3;% kg/m^3
    cv =0.22*u.W/u.KG*rho;
    Q=0; % uniform heat source-- none
    Tampl=100*u.K;
    Tamb=0;
    Tbar =Tampl;%hot face temperature
    num_integ_pts=2; % quadrature
    L=0.6*u.M;% thickness
    W=L/15;% In-plane dimension of the cross-section
    gtolerance=W/1000;
    dt=0.05; % time step
    tend=30; % length of the time interval
    t=0;
    theta = 1.0; % generalized trapezoidal method
    online_graphics= true;% plot the solution as it is computed?
    n=3*5;% needs to be multiple of five
    
    [fens,fes] = H8_block(L,W,W,n,1,1); % Mesh
    
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.conductivity =kappa;
    region.specific_heat =cv;
    region.fes= fes;
    region.integration_rule =gauss_rule(struct('dim',3,'order',2));
    model_data.region{1} =region;
    
    clear essential
    essential.temperature=Tbar;
    essential.node_list = fenode_select(fens,struct('box' , bounding_box ([0,0,0;0,W,W]), 'inflate',gtolerance ));
    model_data.boundary_conditions.essential{1} = essential;

    clear initial_condition
    initial_condition.temperature=Tamb;
    model_data.initial_condition =initial_condition;
    
    model_data.observer = @observer;
    model_data.dt= dt;
    model_data.tend= tend;
    
    Tfifth = [];
    model_data =heat_diffusion_transient(model_data);
    
    function observer (t,  model_data)
        Tn =gather_sysvec(model_data.temp);
        if online_graphics
            plot([0, 0],[0,Tampl],'.'); hold on
            plot(model_data.geom.values(:,1),model_data.temp.values,'rx');
            labels('Distance  through the thickness', 'Temperature [Deg Celsius]');
            figure(gcf);
            title (['Time =' num2str(t)]); hold off; pause(0.1);
        end
        Tfifth = [Tfifth Tn(n/5+1)];
    end
end