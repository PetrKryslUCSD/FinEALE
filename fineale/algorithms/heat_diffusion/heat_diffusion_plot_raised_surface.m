function model_data=heat_diffusion_plot_raised_surface(model_data)
% Plot the temperature as a raised surface above a two-dimensional region.
%
% function model_data=heat_diffusion_plot_raised_surface(model_data)
%
% Arguments:
% model_data=model data as returned  by the solution algorithm
% model_data.postprocessing = optional struct with optional fields
%     gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in model_data.gv
%     z_scale=scale of the Z axis, default is empty; if supplied as a real 
%           number, the data aspect ratio of the plot is changed.
%     camera = camera settings vector, as obtained by camget(gv)
%
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
    
    z_scale = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'z_scale'))
            z_scale = model_data.postprocessing.z_scale;
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
    range = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'range'))
            range = model_data.postprocessing.range;
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
    if (isempty(range))
         range =[min(T),max(T)];;
    end
    dcm=data_colormap(struct('range',range,'colormap',jet));
    
    % Create the color field
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
        map_data(dcm, T)));
    
    % Create the 3-D geometry of the surfaces
    geomT=nodal_field(struct ('name', ['geomT'], ...
        'data',[geom.values, temp.values]));
        
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if (draw_mesh), ec='k';
        else, ec='none';
        end
        draw(region.fes, gv, struct ('x',geomT, 'u',0*geomT,...
            'colorfield',colorfield, 'edgecolor',ec,'shrink',1.0));
        if (draw_mesh)
            draw(region.fes, gv, struct ('x',geom, 'u',0*geom, ...
                'facecolor','none'));
        end
    end
    
    if (~isempty(z_scale))
        
        set(gca,'xlim', [min(geom.values(:,1)),max(geom.values(:,1))]);
        set(gca,'ylim', [min(geom.values(:,2)),max(geom.values(:,2))]);
        set(gca,'zlim', [min([0,dcm.range]),max([0,dcm.range])]);
        set(gca,'DataAspectRatio', [1, 1, 1/z_scale])
        
    end
    
    draw_colorbar(gv,struct('colormap',dcm.colormap,...
                'position',[0.86, 0.1, 0.025, 0.5],...
                'minmax', dcm.range,...
                'label','Temp.', 'fontname', 'Times', 'interpreter', 'latex'));
                
            
    xlabel('X [length]')
    ylabel('Y [length]')
    zlabel('Temperature [degree]')
    
    % Set the camera if supplied
    if (~isempty( camera ))
        camset(gv,camera);
    end
    
    % Return model_data for further use
    model_data.postprocessing.gv=gv;
end