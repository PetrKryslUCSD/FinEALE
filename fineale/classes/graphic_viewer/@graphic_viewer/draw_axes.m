function draw_axes(self, context)
% Draw the basis vectors of the global Cartesian coordinate system.
%
% function draw_axes(self, context)
%
% Draw the basis vector of the global Cartesian coordinate system as
% arrows, red, green, and blue.
% 
% context= optional structure, with fields
%    length = length of the axes' arrows
%    origin = origin of the decoration that represents the coordinate axes
%    length_units= number, representing the length units in which the
%         graphic should be displayed
% 
    if ((nargin >1) && isstruct(context) && isfield(context,'length'))
        a=context.length;
    else
        a=abs(get(self.axes,'XLim'));
        a=a+abs(get(self.axes,'YLim'));
        a=a+abs(get(self.axes,'ZLim'));
        a=sum(a)/20;
        context = struct('color', [0, 0, 0]);
    end
    origin=[0 0 0];
    if ((nargin >1) && isstruct(context) && isfield(context,'origin'))
        origin=context.origin;
    end
    context.color='red';
    draw_arrow(self, origin+[0 0 0], [a  0 0],context);
    draw_text(self, origin+[a  0 0], 'X', context);
    context.color='green';
    draw_arrow(self, origin+[0 0 0], [0  a 0], context);
    draw_text(self, origin+[0 a  0], 'Y', context);
    context.color='blue';
    draw_arrow(self, origin+[0 0 0], [0  0 a], context);
    draw_text(self, origin+[0 0  a], 'Z', context);
    return;
end

