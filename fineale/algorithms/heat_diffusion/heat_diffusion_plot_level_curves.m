function model_data=heat_diffusion_plot_level_curves(model_data)
% Plot the temperature as level curves for a two-dimensional region.
%
% function model_data=heat_diffusion_plot_level_curves(model_data)
%
% Arguments:
% model_data=model data as returned  by the solution algorithm
% model_data.postprocessing= struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in model_data.postprocessing.gv
%      z_scale=scale of the Z axis, default is empty; if supplied as a real 
%           number, the data aspect ratio of the plot is changed.
%      color=  color of the level curve, if not supplied it is computed 
%           from the data
%
% Output
% model_data = structure on input updated with
%  .gv=graphic viewer used to display the data

    
    z_scale = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'z_scale'))
            z_scale = model_data.postprocessing.z_scale;
        end
    end
    color = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'color'))
            color = model_data.postprocessing.color;
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
    dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',jet));
    
    % Create the color field
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
        map_data(dcm, T)));
    
    % Create the 3-D geometry of the surfaces
    geomT=nodal_field(struct ('name', ['geomT'], ...
        'data',[geom.values, temp.values]));
        
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        draw(region.fes, gv, struct ('x',geom, 'u',0*geom, ...
            'facecolor','none'));
        for  isovalue=linspace(min(T),max(T),20)
            if (isempty(color)), c=map_data(dcm, isovalue);
            else, c=color;
            end
            draw_isosurface(region.fes,gv, struct ('x', geomT,'u',0*geomT,'scalarfield',temp,'isovalue',isovalue,'color',c));
        end
    end
    
    if (~isempty(z_scale))
        
        set(gca,'xlim', [min(geom.values(:,1)),max(geom.values(:,1))]);
        set(gca,'ylim', [min(geom.values(:,2)),max(geom.values(:,2))]);
        set(gca,'zlim', dcm.range);
        set(gca,'DataAspectRatio', [1, 1, 1/z_scale])
        
    end
    
    draw_colorbar(gv,struct('colormap',dcm.colormap,...
                'position',[0.86, 0.1, 0.025, 0.5],...
                'minmax', dcm.range,...
                'label','Temp.', 'fontname', 'Times', 'interpreter', 'latex'));
                
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Temperature [degrees C]')
    
    % Return options for further use
    model_data.postprocessing.gv=gv;
end

