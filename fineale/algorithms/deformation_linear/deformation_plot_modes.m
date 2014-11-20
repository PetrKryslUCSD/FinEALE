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
%      save_frame= should we save images for the modes displayed? default false;
%      frame_name= name for the mode images; default 'frame_name';
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
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, u_magnitude.values)));
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
                %     For speed, the deformation is drawn using only the
                %     surface of the solid
                boundaryfes = mesh_boundary (region.femm.fes,[]);
                draw(boundaryfes, gv, struct ('x',geom, 'u',0*model_data.u, ...
                    'facecolor','none'));
                draw(boundaryfes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                    'colorfield',colorfield, 'shrink',1.0));
                %                 The code below would draw the entire mesh, including the interior.
                %                 draw(region.femm.fes, gv, struct ('x',geom, 'u',xscale*model_data.u,...
                %                     'facecolor','none'));
            end
%             draw_text(gv, annloc,['  Eigenvector ' num2str(i) '\newline     ' num2str(sqrt(Omega(ix(i),ix(i)))/2/pi) ' [Hz]'],...
%                 struct('fontsize', 18));
            axis off
            draw_axes (gv,struct('length',axis_length));
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
