function model_data=deformation_plot_steady_vibration(model_data)
% Plot the harmonic force deformation of the structure.
%
% function model_data=deformation_plot_steady_vibration(model_data)
%
% Arguments
% model_data= model data as produced by deformation_linear_statics()
% model_data.postprocessing= optional struct with optional fields
%      gv = graphic viewer; if not supplied, a graphic viewer is created 
%           and returned in options.gv
%      u_scale = deflection scale, default 1.0;
%      frequencylist= default is 1:model_data.neigvs.
%      save_movie= should we save images for the animation frames displayed? default false;
%      movie_name= name for the frame images; default 'movie';
%      draw_mesh= should the mesh be rendered?  Boolean.  Default false.
%      camera  = camera, default is [] which means use the default orientation 
%           of the view;
%      cmap= colormap (default: jet)
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
    camera  = [];
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'camera'))
            camera = model_data.postprocessing.camera;
        end
    end
    frequencylist=1:length(model_data.frequencies);
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'frequencylist'))
            frequencylist = model_data.postprocessing.frequencylist;
        end
    end
    save_movie  =false;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'save_movie'))
            save_movie = model_data.postprocessing.save_movie;
        end
    end
    delay=20;
    movie_name  = 'movie';
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'movie_name'))
            movie_name =  model_data.postprocessing.movie_name;
        end
    end
    cmap = jet;
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
    draw_mesh  = false;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'draw_mesh'))
            draw_mesh = model_data.postprocessing.draw_mesh;
        end
    end
    ncycles  = 1;
    if (isfield(model_data, 'postprocessing'))
        if (isfield(model_data.postprocessing, 'ncycles'))
            ncycles = model_data.postprocessing.ncycles;
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
    range=[inf,-inf]; 
    for i=frequencylist
        % Scatter the eigenvector   
        model_data.u.values=  model_data.u_values{i};
        % Update the range of the amplitudes
        u_magn =magnitude(model_data.u);
        umi=min(u_magn.values);
        uma=max(u_magn.values);
        range=[min([range(1),umi]), max([range(2),uma])]
    end
    max_u_magn = max( range );
    
    for Jj=frequencylist
        frame=1;
        frequency=model_data.frequencies(Jj);% harmonic forcing
        clf; gv=clear (gv,[]); gv=reset (gv,[]);
        if (~isempty(add_to_scene))
            gv=add_to_scene(gv);
        end
        disp(['Frequency ' num2str(Jj) '= ' num2str(frequency) ' Hz']);
        model_data.u.values= model_data.u_values{Jj};
        if (~isempty(map_to_color_fun))
            thecolors=map_to_color_fun(model_data.u);
        else
            % Default: Create the color field from the magnitude of the
            % displacement field
            temp =magnitude(model_data.u); u_magn=temp.values; clear temp
            dcm=data_colormap(struct('range',[min(abs(u_magn)),max(abs(u_magn))],'colormap',cmap));
            thecolors=nodal_field(struct ('name', ['thecolors'], 'data',map_data(dcm, u_magn)));
        end
        % Compute the physical scale to be used for display
        phscale=u_scale*Characteristic_dimension/10*abs(log10(max_u_magn))/abs(log10(norm(u_magn,inf)))*max_u_magn/norm(u_magn,inf);
        ReW=real(model_data.u_values{Jj});
        ImW=imag(model_data.u_values{Jj});
        phases=(0:1:(21*ncycles-1))/21*2*pi;
        for k123=1:length(phases)
            ph=phases(k123);
            model_data.u.values= cos(ph)*ReW-sin(ph)*ImW;
            gv=reset (gv,[]);
            if (~isempty(add_to_scene))
                gv=add_to_scene(gv);
            end
            if (~isempty( camera ))
                camset(gv,camera);
            end
            % Plot the surface for each region
            for i=1:length(model_data.region)
                region =model_data.region{i};
                %     For speed, the deformation is drawn using only the surface of the solid
                boundaryfes = mesh_boundary (region.femm.fes,[]);
                if (draw_mesh )
                    draw(boundaryfes, gv, struct ('x',geom, 'u',0*model_data.u, ...
                        'facecolor','none'));
                end
                if (strcmp(class(thecolors),'nodal_field'))
                    draw(boundaryfes, gv, struct ('x',geom, 'u',phscale*model_data.u,...
                        'colorfield',thecolors, 'edgecolor','none','shrink',1.0));
                else
                draw(boundaryfes, gv, struct ('x',geom, 'u',phscale*model_data.u,...
                        'facecolor',thecolors, 'edgecolor','none','shrink',1.0));
                end
            end
            set_presentation_defaults
            axis off
            if (~isempty(add_decorations))
                gv=add_decorations(gv,frequency,axis_length);
            else
            draw_annotation(gv, [0.01,0.01,0.6,0.1], ...
                ['Freq=' num2str(frequency) ' Hz'], ...
                    struct('color','w','backgroundcolor','w','fontsize',20));
                draw_axes (gv,struct('length',axis_length));
            end
            headlight(gv);
            pause(0.05);
            % Interact with the plot
            interact(gv);
            camera=camget(gv);
            pause(0.05);
            if save_movie % Should we save a movie?
                gif_animation_add_frame(gv.figure,frame,[movie_name '-f' num2str(Jj) '.gif'],40);
                %                 imn =frame_name('.',[movie_name '-f' num2str(Jj) '-'],frame,'png');
                %                 saveas(gv.figure, imn,'png');
            end
            frame=frame+1;
        end
        
        % Interact with the plot
        %         interact(gv);
        
        % Return options for further use
        model_data.postprocessing.gv=gv;
    end
