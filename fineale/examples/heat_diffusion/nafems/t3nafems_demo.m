% T3NAFEMS benchmark problem.
function t3nafems_demo
    kappa=[35.0]; % conductivity matrix
    cm = 440.5;% specific heat per unit mass
    rho=7200;% mass density
    cv =cm* rho;% specific heat per unit volume
    Tampl=100;
    Tamb=0;
    Tbar =@(t)(Tampl*sin(pi*t/40)+ Tamb);%hot face temperature
    L=0.1;% thickness
    dt=2; % time step
    tend= 8*32; % length of the time interval
    online_graphics= ~true;% plot the solution as it is computed?
    n=2*5;% needs to be multiple of five
    
    [fens,fes] = L2_block(L,n,1.0); % Mesh
    
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.conductivity =kappa;
    region.specific_heat =cv;
    region.fes= fes;
    region.integration_rule =trapezoidal_rule(struct('dim',1));
    model_data.region{1} =region;
    
    clear essential
    essential.temperature=Tbar;
    essential.node_list = 1;
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.temperature=Tamb;
    essential.node_list = n+1;
    model_data.boundary_conditions.essential{2} = essential;
    
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
        Tfifth = [Tfifth Tn(n/5)];
        if online_graphics
            plot(((1:length(Tfifth))-1)*dt,Tfifth,'linewidth', 3)
            figure(gcf);
            set(gca,'FontSize', 14); grid on
            xlabel('Time in seconds')
            ylabel('Temperature in degrees Celsius')
            title (['Time =' num2str(t)]); hold off; pause(0.1);
        end
    end
    
    if ~online_graphics
        plot(((1:length(Tfifth))-1)*dt,Tfifth,'linewidth', 3)
        figure(gcf);
        set(gca,'FontSize', 14); grid on
        xlabel('Time in seconds')
        ylabel('Temperature in degrees Celsius')
        title (['Time =' num2str(tend)]);
    end
    
end