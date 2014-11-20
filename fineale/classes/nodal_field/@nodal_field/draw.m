function draw (self, gv, context)
% Produce a graphic representation of the field.
%
% function draw (self, gv, context)
%
%
% Input arguments
% self = self
% gv = graphic viewer
% context = struct
% with mandatory fields
%    x=reference geometry field
%    u=displacement field
%    dof=direction for which to draw the fixed value arrow
%    length=length of the arrows
% with optional fields
%    color = color of the arrows
%    show_free= Boolean, if true show the degrees of freedom which are free 
%
    x = context.x.values; % coordinates of nodes
    u = context.u.values; % Field values at nodes
    direction = zeros(1,length(context.x.values(1, :)));
    direction(context.dof) = context.length;
    show_free=0;
    if (isfield(context,'show_free'))
        show_free=context.show_free;
    end
    if (show_free)
        for i=1:size(self.is_fixed, 1)
            if (~self.is_fixed(i, context.dof))
                draw_arrow(gv, x(i,:)+u(i,:), direction, context);
            end
        end
    else
        for i=1:size(self.is_fixed, 1)
            if (self.is_fixed(i, context.dof))
                draw_arrow(gv, x(i,:)+u(i,:), direction, context);
            end
        end
    end
end