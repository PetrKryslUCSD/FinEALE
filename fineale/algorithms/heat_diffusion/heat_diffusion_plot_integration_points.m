function model_data=heat_diffusion_plot_integration_points(model_data)
% Plot the heat flux at the integration points.
%
% function options=heat_diffusion_plot_integration_points(model_data)
%
% Arguments:
% model_data=model data as returned  by the solution algorithm
% model_data.postprocessing= optional struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in options.gv
%      scale=scale of the integration-point heat flux vectors
%      color=  color of the vectors, if not supplied it is 'red'
%      camera = camera settings vector, as obtained by camget(gv)
%
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
    
    scale  = 1.0;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'scale'))
            scale = model_data.postprocessing.scale;
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
    outputRm  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'outputRm'))
            outputRm = model_data.postprocessing.outputRm;
        end
    end
    
    temp =model_data.temp;
    geom =model_data.geom;
    T=temp.values;
    
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
    
    % Plot the surface for each region
    context=struct ('x',geom,'u',0*geom, ...
            'theta', temp, 'scale', scale,'color','red');
    if (~isempty(outputRm))    
        context.outputRm=outputRm;
    else
        context.outputRm=eye(geom.dim);
    end    
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if (draw_mesh), ec='k';
        else, ec='none';
        end
        draw_integration_points(region.femm, gv, context);    
        if (draw_mesh)
            draw(region.fes, gv, struct ('x',geom, 'u',0*geom, ...
                'facecolor','none'));
        end
    end
    
    xlabel('X ')
    ylabel('Y ')
    
    % Set the camera if supplied
    if (~isempty( camera ))
        camset(gv,camera);
    end
    
    % Return options for further use
    model_data.postprocessing.gv=gv;
end