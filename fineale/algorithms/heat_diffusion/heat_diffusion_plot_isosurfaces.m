function model_data=heat_diffusion_plot_isosurfaces(model_data)
% Plot the temperature as isosurfaces for a three-dimensional region.
%
% function  model_data=heat_diffusion_plot_isosurfaces(model_data)
%
% Arguments:
% model_data=model data as returned  by the solution algorithm
% model_data.postprocessing= optional struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in model_data.postprocessing.gv
%      color=  color of the isosurface, if not supplied it is computed 
%           from the data
%      isovalues= array of values for the iso-surfaces, if not 
%           supplied 20 surfaces are plotted between the minimum and 
%           maximum temperature
%      colormap= colormap (default is jet)
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data

    
    color = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'color'))
            color = model_data.postprocessing.color;
        end
    end
    
    isovalues = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'isovalues'))
            isovalues = model_data.postprocessing.isovalues;
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
    
    % Was the list of iso-values supplied?
      % If not, generate it.
    if (isempty(isovalues))
        isovalues=linspace(min(T),max(T),20);
    end
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        draw(region.fes, gv, struct ('x',geom, 'u',0*geom, ...
            'facecolor','none'));
        for  isovalue=isovalues
            if (isempty(color)), c=map_data(dcm, isovalue);
            else, c=color;
            end
            draw_isosurface(region.fes,gv, struct ('x', geom,'u',0*geom,'scalarfield',temp,'isovalue',isovalue,'color',c));
        end
    end
    
    draw_colorbar(gv,struct('colormap',dcm.colormap,...
                'position',[0.86, 0.1, 0.025, 0.5],...
                'minmax', dcm.range,...
                'label','Temp.', 'fontname', 'Times', 'interpreter', 'latex'));
                
    xlabel('X [length]')
    ylabel('Y [length]')
    zlabel('Z [length]')
    
    % Return options for further use
    model_data.postprocessing.gv=gv;
end

