function model_data=deformation_plot_stress_elementwise(model_data)
% Plot the elementwise stress in the structure.
%
% function model_data=deformation_plot_stress_elementwise(model_data)
%
% Arguments
% model_data= model data as produced by deformation_linear_statics() or by
% a nonlinear solver, for instance deformation_nonlinear_statics().
% model_data.postprocessing = optional struct with optional attributes
%      gv = graphic viewer; if not supplied, a graphic viewer is created
%           and returned in model_data.postprocessing.gv
%      u_scale = deflection scale, default 1.0;
%     output='Cauchy', 'princCauchy', 'pressure', 'vm'
%     stress_range= default is []; For default value the stress range
%           is computed from the stress in the first region.
%     stress_component= default 1 (even for output='pressure');
%     outputRm=  output the stress in this coordinate system
%     stress_units=  physical units for the stress (e.g. u.MEGA*u.PA).
%     camera  = camera, default is [] which means use the default
%           orientation of the view;
%     cmap=colormap, default is jet
%     add_to_scene= function handle, function with signature
%             function gv=add_to_scene(gv);
%          which can be used to add graphics to the viewer (such as spatial cues, or
%          immovable objects)
%     map_to_color_fun= function handle, function with signature
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
%     add_decorations=function handle, function with signature
%             gv=add_decorations(gv,cmap);
%          where gv= graphic viewer, cmap= current colormap.
%          Default: show the color bar and the axes.
%
% Output
% model_data = structure on input is returned updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
%

u_scale = 1.0;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'u_scale'))
        u_scale = model_data.postprocessing.u_scale;
    end
end
stress_range= [];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_range'))
        stress_range = model_data.postprocessing.stress_range;
    end
end
output= 'Cauchy';
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'output'))
        output = model_data.postprocessing.output;
    end
end
stress_component= 1;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_component'))
        stress_component = model_data.postprocessing.stress_component;
    end
end
outputRm= [];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'outputRm'))
        outputRm = model_data.postprocessing.outputRm;
    end
end
stress_units= 1;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_units'))
        stress_units = model_data.postprocessing.stress_units;
    end
end
camera  = [];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'camera'))
        camera = model_data.postprocessing.camera;
    end
end
cmap = jet;
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

% Retrieve displacement fields: either for a linear model (u) or for a nonlinear model
if (isfield(model_data,'u'))
    un1 =model_data.u;    un =0*model_data.u;  dt=[];
elseif (isfield(model_data,'un1'))
    un1 =model_data.un1;    un =model_data.un;  dt =model_data.dt;
end

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
    umax= norm( model_data.u,inf);
    gv=reset (gv,struct('limits', inflate_box(bounding_box(model_data.fens.xyz),u_scale*umax)));
    set_graphics_defaults;
end

% Create the color mapping
if (~isempty(stress_range))
    dcm=data_colormap(struct('range',stress_range,'colormap',cmap));
end

% draw(femm,gv, struct ('x', geom3, 'u', +scale*u3,'colorfield',colorfield, 'edgecolor', 'black','shrink',1.0));
% %     draw(femm,gv, struct ('x', geom3, 'u', +0*u3,'facecolor','none', 'edgecolor','black'));
% draw_colorbar(gv, struct('colormap',cmap,'position',[0.85 0.15 0.05 0.7],...
% 'minmax',nvalsrange,'label',['\sigma_{' num2str(sigj) '}']));
% % view (2)

% Should we create context to pass additional arguments?
context =[];
if (~isempty(outputRm))
    context.outputRm=outputRm;
end


% Plot the surface for each region
minmin=inf; maxmax=-inf;
for i=1:length(model_data.region)
    region =model_data.region{i};
    % Find the finite elements that are exposed on the boundary
    boundaryfes = mesh_boundary (region.femm.fes,[]);
    oon_boundary =zeros(count(model_data.fens),1);
    oon_boundary(boundaryfes.conn(:)) =1;
    for feix=1:count(model_data.region{i}.femm.fes)
        region.femm.fes =subset(model_data.region{i}.femm.fes,feix);
        if (~boundary_only) || (boundary_only)&&(sum(oon_boundary(region.femm.fes.conn(:)))>0)
            % Create the color field
            fld = field_from_integration_points (region.femm, ...
                geom, un1, un, dt, [], output,stress_component, context);
            fld = (1/stress_units)*fld;
            minmin=min([minmin,min(fld.values)]);
            maxmax=max([maxmax,max(fld.values)]);
            if (isempty(stress_range))
                stress_range =[min(fld.values),max(fld.values)];
                dcm=data_colormap(struct('range',stress_range,'colormap',cmap));
            end
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
                map_data(dcm, fld.values)));
            draw(region.femm, gv, struct ('x',geom, 'u',u_scale*un1,...
                'colorfield',colorfield, 'shrink',1.0));
        end
    end
end

% Set the camera if supplied
if (~isempty( camera ))
    camset(gv,camera);
end

% If desired, add additional graphics to the scene
if (~isempty(add_to_scene))
            gv=add_to_scene(gv);
end
      
% If desired, add decorations (annotation, title, ...)
if (~isempty(add_decorations))
    gv=add_decorations(gv,dcm);
else
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    draw_colorbar(gv,struct('colormap',dcm.colormap,...
        'position',[0.81, 0.1, 0.025, 0.5],...
        'minmax', dcm.range,...
        'label',['$\sigma_{' num2str(stress_component) '}$'], 'fontname', 'Times', 'interpreter', 'latex'));
end


% Interact with the plot
interact(gv);

% Return options for further use
model_data.postprocessing.gv=gv;
end
