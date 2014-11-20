function model_data=deformation_plot_integration_points(model_data)
    % Plot the integration point quantities in the structure.
    %
    % function model_data=deformation_plot_integration_points(model_data)
    %
    % Arguments
    % model_data= model data as produced by deformation_linear_statics()
    % model_data.postprocessing = optional struct with optional attributes
    %      gv = graphic viewer; if not supplied, a graphic viewer is created
    %           and returned in model_data.postprocessing.gv
    %      u_scale = deflection scale, default 1.0;
    %      stress_range= default is []; For default value the stress range
    %           is computed from the stress in the first region.
    %      output='Cauchy', 'princCauchy', 'pressure'
    %      stress_component= default 1 (even for output='pressure');
    %      outputRm=  output the stress in this coordinate system
    %      camera  = camera, default is [] which means use the default orientation
    %           of the view;
    %      cmap= colormap (default: jet)
    %
    % Output
    % model_data = structure on input is returned updated with
    % model_data.postprocessing.gv=graphic viewer used to display the data
    %
    
    u_scale = 1.0;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'u_scale'))
            u_scale = model_data.postprocessing.u_scale;
        end
    end
    scale = 1.0;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'scale'))
            scale = model_data.postprocessing.scale;
        end
    end
    stress_range= [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'stress_range'))
            stress_range = model_data.postprocessing.stress_range;
        end
    end
    output= 'Cauchy';
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'output'))
            output = model_data.postprocessing.output;
        end
    end
    stress_component= 1;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'stress_component'))
            stress_component = model_data.postprocessing.stress_component;
        end
    end
    outputRm= [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'outputRm'))
            outputRm = model_data.postprocessing.outputRm;
        end
    end
    draw_mesh = true;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'draw_mesh'))
            draw_mesh = model_data.postprocessing.draw_mesh;
        end
    end
    camera  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'camera'))
            camera = model_data.postprocessing.camera;
        end
    end
    cmap = jet;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'cmap'))
            cmap = model_data.postprocessing.cmap;
        end
    end
     
    u =model_data.u;
    geom =model_data.geom;
    dT=[];
    
    % Initialize the graphic viewer
    gv = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'gv'))
            gv = model_data.postprocessing.gv;
        end
    end
    if (isempty(gv))
        gv=graphic_viewer;
        gv=reset (gv,[]);
        set_graphics_defaults;
    end
    
    % Create the color mapping
    dcm= [];
    if (~isempty(stress_range))
        dcm=data_colormap(struct('range',stress_range,'colormap',cmap));
    end
    
    % Should we create context to pass additional arguments?
    context.x=geom;
    context.u=u;
    context.dT=dT;
    context.quantity = output;
    context.scale=scale;
    context.data_cmap=dcm;
    context.component=stress_component;
    if (~isempty(outputRm))
        context.outputRm=outputRm;
    end
    
    
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if (draw_mesh), ec='k';
        else, ec='none';
        end
        draw_integration_points(region.femm, gv, context);
        if (draw_mesh)
            draw(region.fes, gv, struct ('x',geom, 'u',u_scale*u, ...
                'facecolor','none'));
        end
    end
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    % Set the camera if supplied
    if (~isempty( camera ))
        camset(gv,camera);
    end
    
    if (~isempty(dcm))
        draw_colorbar(gv,struct('colormap',dcm.colormap,...
            'position',[0.86, 0.1, 0.025, 0.5],...
            'minmax', dcm.range,...
            'label',['$\sigma_{' num2str(stress_component) '}$'], 'fontname', 'Times', 'interpreter', 'latex'));
    end
    
    % Interact with the plot
    interact(gv);
    
    % Return options for further use
    model_data.postprocessing.gv=gv;
end
