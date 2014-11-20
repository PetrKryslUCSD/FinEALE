function model_data=heat_diffusion_plot_temperature(model_data)
% Plot the temperature as a color-coded surface.
%
% function model_data=heat_diffusion_plot_temperature(model_data)
%
% Arguments:
% model_data=model data as returned  by the solution algorithm
% model_data.postprocessing= optional struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in model_data.postprocessing.gv
%      camera = camera settings vector, as obtained by camget(gv)
%      colormap= colormap (default is jet)
%      draw_mesh= should the mesh be rendered?  Boolean.  Default true.
%
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
    
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
    colormap = jet;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'colormap'))
            colormap = model_data.postprocessing.colormap;
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
    
    % Create the color mapping
    dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',colormap));
    
    % Create the color field
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
        map_data(dcm, T)));
    
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if (draw_mesh), ec='k';
        else, ec='none';
        end
        draw(region.fes, gv, struct ('x',geom, 'u',0*geom,...
            'colorfield',colorfield, 'edgecolor',ec,'shrink',1.0));
    end
    
    draw_colorbar(gv,struct('colormap',dcm.colormap,...
                'position',[0.86, 0.1, 0.025, 0.5],...
                'minmax', dcm.range,...
                'label','Temp.', 'fontname', 'Times', 'interpreter', 'latex'));
                
            
    xlabel('X [length]')
    ylabel('Y [length]')
    zlabel('Z [length]')
    
    % Set the camera if supplied
    if (~isempty( camera ))
        camset(gv,camera);
    end
    
    % Return options for further use
    model_data.postprocessing.gv=gv;
end