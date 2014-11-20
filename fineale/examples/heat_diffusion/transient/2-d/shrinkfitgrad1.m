% Temperature equalization  in shrink-fitting. Uniform mesh.
function shrinkfitgrad1
kappa_steel=[44 0; 0 44]; % N/K/s
kappa_tungsten=[163 0; 0 163]; % N/K/s
rho_steel = 7500*1e-9;% kg/mm^3
rho_tungsten = 19000*1e-9;% kg/mm^3
cv_steel =470*1e3*rho_steel;
cv_tungsten =134*1e3*rho_tungsten;
h=50*1e-3;
Ta=17;
Ts=84;
Tt=-10;
dt=2.5; % time step, try also 5.0, 2.0
theta =1.0; % generalized trapezoidal method parameter
tend= 100; % length of the time interval
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    'curve 1 line 0 0 50 0',...
    'curve 2 arc 50 0 80 0 center 65 -0.001 ',...
    'curve 3 line 80 0 110 0',...
    'curve 4 line 110 0 110 50',...
    'curve 5 line 110 50 65 50 ',...
    'curve 6 arc 65 50 65 70 center 65.001 60  ',...
    'curve 7 line 65 70 110 70',...
    'curve 8 line 110 70 110 85',...
    'curve 9 arc 110 85 65 120 center 110 120 ',...
    'curve 10 line 65 120 0 120',...
    'curve 11 line 0 120 0 85',...
    'curve 12 arc 0 85 0 35 center -0.001 60 rev',...
    'curve 13 line 0 35 0 0',...
    'curve 14 line 110, 50, 160, 50',...
    'curve 15 line 160, 50, 160, 70',...
    'curve 16 line 160, 70, 110, 70',...
    ['subregion 1  property 1 boundary '...
    ' 1 2 3 4 5 6 7 8 9 10 11 12 13'],...
    ['subregion 2  property 2 boundary '...
    ' -5 -6 -7 14 15 16'],...
    ['m-ctl-point constant 7.']
    }, 1.0);


    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.conductivity =kappa_steel;
    region.specific_heat =cv_steel;
    region.fes= subset(fes,groups{1});
    region.integration_rule =tri_rule(struct('npts',1));
    model_data.region{1} =region;
    
    clear region
    region.conductivity =kappa_tungsten;
    region.specific_heat =cv_tungsten;
    region.fes= subset(fes,groups{2});
    region.integration_rule =tri_rule(struct('npts',1));
    model_data.region{2} =region;
    
    clear convection
    convection.ambient_temperature=Ta;
    convection.surface_transfer_coefficient  = h;
    convection.fes = subset(edge_fes,[edge_groups{[(1:4) (8:16)]}]);;
    convection.integration_rule= gauss_rule(struct('dim',1,'order',1));
    model_data.boundary_conditions.convection{1} = convection;

    function Ti = initial_temperature (xx)
    Ti =zeros(size(xx,1),1)+Ts;
    conn = connected_nodes(model_data.region{2}.fes);
    Ti(conn) = Tt;
    end
    clear initial_condition
    initial_condition.temperature=@initial_temperature;
    model_data.initial_condition =initial_condition;
    
    model_data.observer = @observer;
    model_data.dt= dt;
    model_data.tend= tend;
    model_data.theta=  theta;
    
    gv=graphic_viewer;
    set_graphics_defaults(gcf)
    axis equal vis3d
    minmaxT = [];
    model_data =heat_diffusion_transient(model_data);
    
    
    figure
    plot((0:1:size(minmaxT, 1)-1)*dt,minmaxT,'linewidth', 3)
    set(gca,'FontSize', 14); grid on
    xlabel('Time in seconds')
    ylabel('Temperature in degrees Celsius')
    
    function observer (t,  model_data)
        gv=reset (gv,[]);
        T=model_data.temp.values;
        minmaxT= [minmaxT; min(T),max(T)]
        model_data.postprocessing.scale = 0.05;
        model_data=heat_diffusion_plot_integration_points(model_data);
        view(2);
        figure(gcf);
        title (['Time =' num2str(t)]); hold off; pause(1);
        %                 saveas(gcf, ['shrinkfit-' num2str(t) '.png'], 'png');
    end
    
    end