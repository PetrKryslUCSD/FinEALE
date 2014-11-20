% Transient simulation to demonstrate the smearing out of the initial temperature 
function hill_decay
    kappa=[0.2 0; 0 0.2]; % conductivity matrix
    cv=1.0; % specific heat per unit volume
    theta = 1/2; % generalized trapezoidal method parameter
    dt=2; % time step
    tend= 20*dt; % length of the time interval
    T0=@(xy)500./(xy(:,1).^2+xy(:,2).^2+5);% Initial distribution of temperature
    
    mesh_size =  2.;
    [fens,fes] = targe2_mesher({...
        ['curve 1 line -30 -20 30 -20'],...
        ['curve 2 line 30 -20 30 20'],...
        ['curve 3 line 30 20 -30 20'],...
        ['curve 4 line -30 20 -30 -20'],...
        'subregion 1  property 1 boundary 1 2 3 4 ',...
        ['m-ctl-point constant '   num2str(mesh_size)]
        }, 1.0);
    
    tolerance =mesh_size/100;
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.conductivity =kappa;
    region.specific_heat =cv;
    region.fes= fes;
    region.integration_rule =tri_rule(struct('npts',1));
    model_data.region{1} =region;
    
    clear initial_condition
    initial_condition.temperature=T0;
    model_data.initial_condition =initial_condition;
    
    model_data.observer = @observer;
    model_data.dt= dt;
    model_data.tend= tend;
    
    gv=graphic_viewer;
    set_graphics_defaults(gcf)
    axis equal vis3d
    
    model_data =heat_diffusion_transient(model_data);
    
    function observer (t,  model_data)
        gv=reset (gv,struct('limits',[-30,30,-20,20,0,100]));
        T=model_data.temp.values;
        geomT=nodal_field(struct ('name', ['geomT'], ...
            'data',[model_data.geom.values, model_data.temp.values]));
        draw(model_data.region{1}.fes, gv, struct ('x',geomT, 'u',0*geomT,...
            'facecolor','y'));
        labels X Y T ;
        camset(gv,[-146.1721 -410.9993  485.8627  -13.8460   -3.7415   42.4319         0         0 1.0000    8.7074]);
        figure(gcf);
        title (['Time =' num2str(t)]); hold off; pause(0.1);
        %                 saveas(gcf, ['shrinkfit-' num2str(t) '.png'], 'png');
    end
    
end