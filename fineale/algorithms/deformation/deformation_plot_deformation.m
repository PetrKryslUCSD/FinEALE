function model_data=deformation_plot_deformation(model_data)
% Plot the deformation of the structure.
%
% function model_data=deformation_plot_deformation(model_data)
%
% Arguments
% model_data= model data as produced by deformation_linear_statics() or by
% a nonlinear solver, for instance deformation_nonlinear_statics().
% model_data.postprocessing = optional struct with optional fields
%     u_scale = deflection scale, default 1.0;
%     quantity='U' (magnitude, default), or 'Un'  (displacement component n)
%     camera  = camera, default is [] which means use the default 
%           orientation of the view;
%     draw_mesh= should the mesh be rendered?  Boolean.  Default true.
%     cmap= colormap (default: jet)
%     boundary_only= should we plot only the boundary?  Boolean.  Default true.
%     edgecolor= color used for the mesh edges
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
%    
    
    u_scale = 1.0;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'u_scale'))
            u_scale = model_data.postprocessing.u_scale;
        end
    end
    quantity = 'U';
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'quantity'))
            quantity = model_data.postprocessing.quantity;
        end
    end
    camera  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'camera'))
            camera = model_data.postprocessing.camera;
        end
    end
    draw_mesh  = true;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'draw_mesh'))
            draw_mesh = model_data.postprocessing.draw_mesh;
        end
    end
    cmap = jet(64);
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'cmap'))
            cmap = model_data.postprocessing.cmap;
        end
    end
    boundary_only= true;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'boundary_only'))
            boundary_only = model_data.postprocessing.boundary_only;
        end
    end
    edgecolor= [1,1,1]*0.;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'edgecolor'))
            edgecolor = model_data.postprocessing.edgecolor;
        end
    end

    colorbar_context.fontname='Times';;
    colorbar_context.interpreter='latex';;
    colorbar_context.position=[0.81, 0.1, 0.025, 0.5];;
    colorbar_context.label=['$u_{mag}$'];;
    colorbar_context.fontname= 'Times';;
    colorbar_context.interpreter= 'latex';
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'colorbar_context'))
            colorbar_context = model_data.postprocessing.colorbar_context;
        end
    end
    
    % Retrieve displacement fields: either for a linear model (u) or for a nonlinear model
    if (isfield(model_data,'u'))
        un1 =model_data.u;    un =0*model_data.u;  dt=[];
    elseif (isfield(model_data,'un1'))
        un1 =model_data.un1;    un =model_data.un;  dt =model_data.dt;
    end

    geom =model_data.geom;
    u_magnitude =magnitude(un1);
    
    % Initialize the graphic viewer
    gv = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'gv'))
            gv = model_data.postprocessing.gv;
        end
    end
    if (isempty(gv))
        gv=graphic_viewer;
        umax= norm( model_data.u,inf);
        gv=reset (gv,struct('limits', inflate_box(bounding_box(model_data.fens.xyz),u_scale*umax)));
        set_graphics_defaults;
    end
    
    
    % Create the color field
    if (strcmp( quantity, 'U1'))
        % Create the color mapping
        dcm=data_colormap(struct('range',[min(un1.values(:,1)),max(un1.values(:,1))],'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, un1.values(:,1))));
    elseif (strcmp( quantity, 'U2'))
        % Create the color mapping
        dcm=data_colormap(struct('range',[min(un1.values(:,2)),max(un1.values(:,2))],'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, un1.values(:,2))));
    elseif (strcmp( quantity, 'U3'))
        % Create the color mapping
        dcm=data_colormap(struct('range',[min(un1.values(:,3)),max(un1.values(:,3))],'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, un1.values(:,3))));
    else
        % Create the color mapping
        dcm=data_colormap(struct('range',[0,max(u_magnitude.values)],'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, u_magnitude.values)));
    end
    
    % Plot the surface for each region
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if boundary_only
            if (draw_mesh)
                draw(mesh_boundary(region.fes,[]), gv, struct ('x',geom, 'u',0*un1, ...
                    'facecolor','none','edgecolor',edgecolor));
            end
            draw(mesh_boundary(region.fes,[]), gv, struct ('x',geom, 'u',u_scale*un1,...
                'colorfield',colorfield, 'shrink',1.0));
        else
            if (draw_mesh)
                draw(region.fes, gv, struct ('x',geom, 'u',0*un1, ...
                    'facecolor','none','edgecolor',edgecolor));
            end
            draw(region.fes, gv, struct ('x',geom, 'u',u_scale*un1,...
                'colorfield',colorfield,'edgecolor',edgecolor, 'shrink',1.0));
        end
    end
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    % Set the camera if supplied
    if (~isempty( camera ))
        camset(gv,camera);
    end
    
    colorbar_context.colormap=dcm.colormap;
    colorbar_context.minmax=dcm.range;;
    draw_colorbar(gv,colorbar_context);
 
    % Interact with the plot
    interact(gv);
    
    % Return model_data for further use
    model_data.postprocessing.gv=gv;
end
