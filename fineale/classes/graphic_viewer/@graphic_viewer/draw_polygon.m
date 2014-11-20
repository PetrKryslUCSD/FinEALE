function draw_polygon(self, x, face, context)
% Draw polygon.
%
% function draw_polygon(self, x, face, context)
%
%
% required arguments:
% x=array of nvert number of rows and 3 columns; the coordinates
% face=array of M rows and N columns; N=the numbers of vertices
%        of an N-sided polygonal patch.  The indices point into
%        the x array.
% context = struct with following optional fields
%    facecolor = color if polygons should be drawn in solid color,
%    colors = array of nvert number of rows and 3 columns; per vertex colors
%    fill=string: either 'none' or 'interp' or 'flat'
%    length_units= number, representing the length units in which the
%         graphic should be displayed
%
    edgecolor='black';
    if isfield(context,'edgecolor')
        edgecolor =context.edgecolor;
    end
    if isfield(context,'length_units')
        x=x/context.length_units;
    end
    if isfield(context,'colors')
        h=patch('Faces', face, ...
            'Vertices', x,'BackFaceLighting','lit',...
            'AmbientStrength',0.75,...
            'FaceVertexCData', context.colors,...
            'FaceColor', 'interp','EdgeColor',edgecolor);%'NormalMode','manual',
    else
        facecolor='red';
        if isfield(context,'facecolor')
            facecolor =context.facecolor;
        end
        h=patch('Faces', face, ...
            'Vertices', x,'BackFaceLighting','lit',...
            'AmbientStrength',0.75,...
            'FaceColor', facecolor,'EdgeColor',edgecolor);%'NormalMode','manual',
    end
    if isfield(context,'facealpha')
        set (h,'FaceAlpha', context.facealpha);
    end
    if isfield(context,'clipping')
        set (h,'clipping', context.clipping);
    end
end
% for j=1:size(face,1)
%             h=patch('Faces', face(j,:), ...
%                 'Vertices', x,'NormalMode','manual','BackFaceLighting','lit',...
%                 'AmbientStrength',0.75,...
%                 'FaceColor', facecolor,'EdgeColor',edgecolor);
%             if isfield(context,'facealpha')
%                 set (h,'FaceAlpha', context.facealpha);
%             end
%         end