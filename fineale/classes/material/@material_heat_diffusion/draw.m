function draw (self, gv, context)
% Produce a graphic representation of data at the integration point.
%
% function draw (self, gv, context)
%
%
% Input arguments
% self = self
% gv = graphic viewer
% context = struct
%  with mandatory fields
%    xyz         - location at which the draw
%    ms          - material state
%    update_context     - context to be passed to the update () function.
%  with optional fields
%    quantity  - string: name of the quantity to display
%    component - index that says which component of the `quantity' quantity;
%                remember: the component points into the output vector
%    data_cmap - data color map to use to color the quantity
%    Rm - transformation 3x3 matrix, with basis vectors in columns, that
%         describes the local Cartesian system for the material point
%
    quantity='flux';
    if (isfield(context,'quantity'))
        quantity= context.quantity;
    end
    scale = 1.0;
    if (isfield(context,'scale'))
        scale= context.scale;
    end
    switch (quantity)
        case 'flux'
            flux = update (self, context.ms, context.update_context);
            if ( isfield(  context,'Rm' ))
                if (~isempty(context.Rm))
                    flux=context.Rm*flux;
                end
            end
            component =1;
            if isfield(context,'component')
                component = context.component;
            end
            color =[1 1 0];
            if isfield(context,'data_cmap')
                color = map_data(context.data_cmap, flux(component));
            end
            context.facecolor= color;
            draw_arrow(gv, context.xyz, scale*flux, context);
        otherwise
            %  do nothing;
    end
    return;
end
