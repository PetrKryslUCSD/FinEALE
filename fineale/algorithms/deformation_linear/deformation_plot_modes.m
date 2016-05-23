function model_data=deformation_plot_modes(model_data)
% Plot the modal deformation of the structure.
%
% function model_data=deformation_plot_modes(model_data)
%
% Arguments
% model_data= model data as produced by deformation_linear_statics()
% model_data.postprocessing= optional struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in options.gv
%      u_scale = deflection scale, default 1.0;
%      modelist= default is 1:model_data.neigvs.
%      draw_mesh= should the mesh be rendered?  Boolean.  Default true.
%      save_frame= should we save images for the modes displayed? default false;
%      frame_name= name for the mode images; default 'frame_name';
%      camera  = camera, default is [] which means use the default orientation 
%           of the view;
%      cmap= colormap (default: jet)
%      animate= should the mode shape be animated? true or false (default)
%      add_to_scene= function handle, function with signature
%             function gv=add_to_scene(gv);
%          which can be used to add graphics to the viewer (such as spatial cues, or
%          immovable objects)
%      map_to_color_fun= function handle, function with signature
%             function v=fun(fld, cmap)
%          where fld= displacement field, cmap=colormap, and the 
%          output v= nodal_field color field (field with three colors per node)
%          or a color  specification (for instance 'y' or [0.8, 0.4, 0.3]).
%          This is optional: default is 
%             temp =magnitude(model_data.u); u_magn=temp.values; clear temp
%             dcm=data_colormap(struct('range',[min(u_magn),max(u_magn)],'colormap',cmap));
%             thecolors=nodal_field(struct ('name', ['thecolors'], 'data',map_data(dcm, u_magn)));
%          The magnitude should be  sqrt(sum(fld.values.*conj(fld.values),2)) for
%          complex-valued fields.
%      add_decorations=function handle, function with signature
%             gv=add_decorations(gv,frequency,axis_length);
%          where gv= graphic viewer, frequency= current frequency, axis_length=
%          an appropriate length for the axes.  Default: display the
%          frequency information with an annotation.
 %
% Output
% model_data = structure on input updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
  
    u_scale = 1.0;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'u_scale'))
            u_scale = model_data.postprocessing.u_scale;
        end
    end
    animate = false;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'animate'))
            animate = model_data.postprocessing.animate;
        end
    end
    camera  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'camera'))
            camera = model_data.postprocessing.camera;
        end
    end
    modelist=1:model_data.neigvs;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'modelist'))
            modelist = model_data.postprocessing.modelist;
        end
    end
    draw_mesh  = true;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'draw_mesh'))
            draw_mesh = model_data.postprocessing.draw_mesh;
        end
    end
    save_frame  = false;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'save_frame'))
            save_frame = model_data.postprocessing.save_frame;
        end
    end
    frame_name  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'frame_name'))
            frame_name =  model_data.postprocessing.frame_name;
        end
    end
    cmap = parula;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'cmap'))
            cmap = model_data.postprocessing.cmap;
        end
    end
    map_to_color_fun = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'map_to_color_fun'))
            map_to_color_fun = model_data.postprocessing.map_to_color_fun;
        end
    end
    add_to_scene  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'add_to_scene'))
            add_to_scene = model_data.postprocessing.add_to_scene;
        end
    end
    add_decorations  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'add_decorations'))
            add_decorations = model_data.postprocessing.add_decorations;
        end
    end
    
    b=update_box([],model_data.fens.xyz);
    axis_length=mean([diff(b(1:2)),diff(b(3:4)),diff(b(5:6))])/4;
    
    geom =model_data.geom;
    
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
    
    % Find the bounding box of the structure and deduce "characteristic" dimension
    Characteristic_dimension  =mean(diff([min(geom.values);max(geom.values)]));
    
    % Create the color mapping
    range=[0,0];
    for i=modelist
        % Scatter the eigenvector   
        model_data.u= scatter_sysvec(model_data.u, model_data.W(:,i));
        % Update the range of the amplitudes
        u_magnitude =magnitude(model_data.u,2);
        range=[min([range(1),min(abs(u_magnitude.values))]),...
            max([range(2),max(abs(u_magnitude.values))])];
    end
    dcm=data_colormap(struct('range',range,'colormap',cmap));
    
    for Jj=modelist
        clf; gv=clear (gv,[]); gv=reset (gv,[]);
        disp(['Mode ' num2str(Jj) ', frequency ' num2str(model_data.Omega(Jj)/2/pi)]);
        model_data.u= scatter_sysvec(model_data.u, model_data.W(:,Jj));
        u_magnitude =magnitude(model_data.u,2);
        % Create the color field
        if (~isempty(map_to_color_fun))
            thecolors=map_to_color_fun(model_data.u);
        else
            % Default: Create the color field from the magnitude of the
            % displacement field
            temp =magnitude(model_data.u); u_magn=temp.values; clear temp
            dcm=data_colormap(struct('range',[min(abs(u_magn)),max(abs(u_magn))],'colormap',cmap));
            thecolors=nodal_field(struct ('name', ['thecolors'], 'data',map_data(dcm, u_magn)));
        end
        if animate
            xscales=u_scale*Characteristic_dimension/10/max( range )*sin((0:1:42)/21*2*pi);
        else
            xscales=u_scale*Characteristic_dimension/10/max( range );
        end
        for xscale=xscales
            gv=reset (gv,[]);
            if (~isempty( camera ))
                camset(gv,camera);
            end
            % Plot the surface for each region
            for i=1:length(model_data.region)
                region =model_data.region{i};
                if (~isempty(add_to_scene))
                    gv=add_to_scene(gv);
                end
                %     For speed, the deformation is drawn using only the
                %     surface of the solid
                boundaryfes = mesh_boundary (region.femm.fes,[]);
                if (draw_mesh)
                    draw(boundaryfes, gv, struct ('x',geom, 'u',0*model_data.u, ...
                        'facecolor','none'));
                end
                if (strcmp(class(thecolors),'nodal_field'))
                draw(boundaryfes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                        'colorfield',thecolors, 'shrink',1.0));
                else
                    draw(boundaryfes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                        'facecolor',thecolors, 'shrink',1.0));
                end
                %                 The code below would draw the entire mesh, including the interior.
                %                 draw(region.femm.fes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                %                     'facecolor','none'));
            end
%             draw_text(gv, annloc,['  Eigenvector ' num2str(i) '\newline     ' num2str(sqrt(Omega(ix(i),ix(i)))/2/pi) ' [Hz]'],...
%                 struct('fontsize', 18));
            axis off
            if (~isempty(add_decorations))
                gv=add_decorations(gv,frequency,axis_length);
            else
                draw_annotation(gv, [0.01,0.01,0.98,0.08], ...
                    ['Mode ' num2str(Jj) ': f=' num2str(model_data.Omega(Jj)/2/pi) ' Hz'], ...
                    struct('color','k','edgecolor','w','backgroundcolor','w','fontsize',18));
            draw_axes (gv,struct('length',axis_length));
            end
            pause(0.005);
        end
        
     
        
        % Set the camera if supplied
        if (~isempty( camera ))
            camset(gv,camera);
        end
        
        if ( save_frame )
            if (isempty(frame_name))
                frame_name ='frame';
            end
            imn =[frame_name '-m' num2str(Jj)];
            draw_annotation(gv, [0.01,0.01,0.8,0.1], ...
            ['Mode ' num2str(Jj) ': Frequency =' num2str(model_data.Omega(Jj)/2/pi) ' Hz'], ...
            struct('color','r','fontsize',20));
            saveas(gcf, imn, 'png');
        end
        
        % Interact with the plot
        interact(gv);
        
        % Return options for further use
        model_data.postprocessing.gv=gv;
    end
