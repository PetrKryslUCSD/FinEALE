function model_data=deformation_plot_stress(model_data)
% Plot the stress in the structure.
%
% function model_data=deformation_plot_stress(model_data)
%
% Arguments
% model_data= model data as produced by deformation_linear_statics()
% model_data.postprocessing = optional struct with optional fields
%     observer = handle of an observer function (optional)
%           The observer function has a signature
%                     function output(i, stressf, model_data)
%           where i= number of the region, stressf=the stress nodal field
%     u_scale = deflection scale, default 1.0;
%     output='Cauchy', 'princCauchy', 'pressure', 'vm'
%     stress_range= default is []; For default value the stress range
%           is computed from the stress in the first region.
%     stress_component= default 1 (even for output='pressure');
%     outputRm=  output the stress in this coordinate system
%     stress_units=  physical units for the stress (e.g. u.MEGA*u.PA).
%     use_spr= should the super convergent patch recovery be used
%           to compute the stress? default false;
%     camera  = camera, default is [] which means use the default
%           orientation of the view;
%     cmap=colormap, default is jet
%
% Output
% model_data  = structure on input is returned updated with
% model_data.postprocessing.gv=graphic viewer used to display the data
%

u_scale = 1.0;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'u_scale'))
        u_scale = model_data.postprocessing.u_scale;
    end
end
output= 'Cauchy';
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'output'))
        output = model_data.postprocessing.output;
    end
end
stress_range= [];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_range'))
        stress_range = model_data.postprocessing.stress_range;
    end
end
stress_component= 1;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_component'))
        stress_component = model_data.postprocessing.stress_component;
    end
end
stress_units= 1;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'stress_units'))
        stress_units = model_data.postprocessing.stress_units;
    end
end
outputRm= [];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'outputRm'))
        outputRm = model_data.postprocessing.outputRm;
    end
end
use_spr= 0;
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'use_spr'))
        use_spr = model_data.postprocessing.use_spr;
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
observer =[];
if (isfield(model_data, 'postprocessing'))
    if (isfield(model_data.postprocessing, 'observer'))
        observer  =model_data.postprocessing.observer;;
    end
end

u =model_data.u;
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

% Create the color mapping
if (~isempty(stress_range))
    dcm=data_colormap(struct('range',stress_range,'colormap',cmap));
end

% Should we create context to pass additional arguments?
context =[];
if (~isempty(outputRm))
    context.outputRm=outputRm;
end

% Plot the surface for each region
for i=1:length(model_data.region)
    region =model_data.region{i};
    % Create the color field
    if (use_spr)
        fld = field_from_integration_points_spr (region.femm, ...
            geom, u, [], output,stress_component,context);
    else
        fld = field_from_integration_points (region.femm, ...
            geom, u, [], output,stress_component, context);
    end
    fld = (1/stress_units)*fld;
    if ~isempty(observer)% report the progress
        observer (i,fld,model_data);
    end
    if (isempty(stress_range))
        stress_range =[min(fld.values),max(fld.values)];
        dcm=data_colormap(struct('range',stress_range,'colormap',cmap));
    end
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
        map_data(dcm, fld.values)));
    if (region.femm.fes.dim>2)
        boundaryfes = mesh_boundary (region.femm.fes,[]);
        draw(boundaryfes, gv, struct ('x',geom, 'u',u_scale*u,...
            'colorfield',colorfield, 'shrink',1.0));
    else
        draw(region.femm.fes, gv, struct ('x',geom, 'u',u_scale*u,...
            'colorfield',colorfield, 'shrink',1.0));
    end
end

xlabel('X')
ylabel('Y')
zlabel('Z')

% Set the camera if supplied
if (~isempty( camera ))
    camset(gv,camera);
end

draw_colorbar(gv,struct('colormap',dcm.colormap,...
    'position',[0.81, 0.1, 0.025, 0.5],...
    'minmax', dcm.range,...
    'label',['$\sigma_{' num2str(stress_component) '}$'], 'fontname', 'Times', 'interpreter', 'latex'));


% Interact with the plot
interact(gv);

% Return model_data for further use
model_data.postprocessing.gv=gv;
end
