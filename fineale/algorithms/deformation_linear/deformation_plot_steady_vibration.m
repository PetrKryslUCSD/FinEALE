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
        u_magnitude =magnitude(model_data.u,2);
        umi=min(abs(u_magnitude.values));
        uma=max(abs(u_magnitude.values));
        range=[min([range(1),umi]), max([range(2),uma])]
    end
    max_u_magn = max( range );
    
    for Jj=frequencylist
        frame=1;
        frequency=model_data.frequencies(Jj);% harmonic forcing
        clf; gv=clear (gv,[]); gv=reset (gv,[]);
        disp(['Frequency ' num2str(Jj) '= ' num2str(frequency) ' Hz']);
        model_data.u.values= model_data.u_values{Jj};
        u_magnitude =magnitude(model_data.u,2);
        u_magn = max(abs(u_magnitude.values));
        dcm=data_colormap(struct('range',[min(abs(u_magnitude.values)),max(abs(u_magnitude.values))],'colormap',cmap));
        % Create the color field
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, u_magnitude.values)));
        phscale=u_scale*Characteristic_dimension/10*abs(log10(max_u_magn))/abs(log10(u_magn))*max_u_magn/u_magn;
        ReW=real(model_data.u_values{Jj});
        ImW=imag(model_data.u_values{Jj});
        phases=(0:1:(21*ncycles-1))/21*2*pi;
        for k123=1:length(phases)
            ph=phases(k123);
            model_data.u.values= cos(ph)*ReW-sin(ph)*ImW;
            gv=reset (gv,[]);
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
                draw(boundaryfes, gv, struct ('x',geom, 'u',phscale*model_data.u,...
                    'colorfield',colorfield, 'edgecolor','none','shrink',1.0));
                %                 The code below would draw the entire mesh, including the interior.
                %                 draw(region.femm.fes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                %                     'facecolor','none'));
            end
            set_presentation_defaults
            axis off
            draw_annotation(gv, [0.01,0.01,0.6,0.1], ...
                ['Freq=' num2str(frequency) ' Hz'], ...
                struct('color','k','backgroundcolor','w','fontsize',20));
            %             draw_axes (gv,struct('length',axis_length));
            headlight(gv);
            pause(0.05);
            % Interact with the plot
            interact(gv);
            camera=camget(gv);
            if save_movie % Should we save a movie?
                imn =frame_name('.',[movie_name '-f' num2str(Jj) '-'],frame,'png');
                saveas(gv.figure, imn,'png');
            end
            frame=frame+1;
        end
        
        % Interact with the plot
        %         interact(gv);
        
        % Return options for further use
        model_data.postprocessing.gv=gv;
    end
