function draw_planes(self, context)
% Draw the coordinate planes of the global Cartesian coordinate system.
%
% function draw_planes(self, context)
%
% Draw the basis vector of the global Cartesian coordinate system as
% arrows, red, green, and blue.
% 
% context= optional structure, with fields
%    length = length of the side of the squares that represent the planes
%    facealpha= transparency
%    origin = origin of the decoration that represents the coordinate axes
%         and others to be interpreted by draw_text()
%    length_units= number, representing the length units in which the
%         graphic should be displayed
    if ((nargin >1) && isstruct(context) && isfield(context,'length'))
        a=context.length;
    else
        a=abs(get(self.axes,'XLim'));
        a=a+abs(get(self.axes,'YLim'));
        a=a+abs(get(self.axes,'ZLim'));
        a=sum(a)/20;
        context = struct('color', [0, 0, 0]);
    end
    facealpha=0.3;
    if (isfield(context,'facealpha'))
        facealpha= context.facealpha;
    end
    Origin=[0 0 0];
    if (isfield(context,'origin'))
        Origin= context.origin;
    end
    context.facecolor='red';
    context.edgecolor='none';
    context.facealpha=facealpha;
    draw_polygon(self, ones(4,1)*Origin+[[0 0 0]; [0 a 0]; [0 a a]; [0 0 a]], [1,2,3,4], context);
    draw_text(self, Origin+[a  0 0], 'X', context);
    context.facecolor='green';
    context.edgecolor='none';
    context.facealpha=facealpha;
    draw_polygon(self, ones(4,1)*Origin+[[0 0 0]; [a  0 0]; [a  0 a]; [0 0 a]], [1,2,3,4], context);
    draw_text(self, Origin+[0 a  0], 'Y', context);
    context.facecolor='blue';
    context.edgecolor='none';
    context.facealpha=facealpha;
    draw_polygon(self, ones(4,1)*Origin+[[0 0 0]; [a  0 0]; [a  a 0]; [0 a 0]], [1,2,3,4], context);
    draw_text(self, Origin+[0 0  a], 'Z', context);
    return;
end

