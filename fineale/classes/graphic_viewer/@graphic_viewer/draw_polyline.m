function draw_polyline(self, x, edge, context)
% Draw polyline.
%
% function draw_polyline(self, x, edge, context)
%
% required arguments:
% x=array of nvert number of rows and 3 columns; the coordinates
% edge=array of N rows and M columns; N=the number of polylines, 
%        M =the numbers of vertices
%        of an N-sided polyLine.  The indices point into
%        the x array.
% context = struct with following optional fields
%    edgecolor = color if polyline should be drawn in solid color,
%    colors = array of nvert number of rows and 3 columns; per vertex colors
%    length_units= number, representing the length units in which the
%         graphic should be displayed
    edgecolor='black';
    if isfield(context,'edgecolor')
        edgecolor =context.edgecolor;
    end
    linewidth=1;
    if isfield(context,'linewidth')
        linewidth =context.linewidth;
    end
    if isfield(context,'length_units')
        x=x/context.length_units;
    end
    if (~strcmp(edgecolor,'none'))
        if (size(x,2)==3)
            X =  reshape(x(edge,1),size(edge,1),size(edge,2));
            Y =  reshape(x(edge,2),size(edge,1),size(edge,2));
            Z =  reshape(x(edge,3),size(edge,1),size(edge,2));
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            %         line(x(edge,1),x(edge,2),x(edge,3), 'Color',edgecolor);
        elseif (size(x,2)==2) 
            X =  reshape(x(edge,1),size(edge,1),size(edge,2));
            Y =  reshape(x(edge,2),size(edge,1),size(edge,2));
            line(X',Y', 'Color',edgecolor,'linewidth',linewidth);
        else
            X =  reshape(x(edge,1),size(edge,1),size(edge,2));
            line(X',0*X', 'Color',edgecolor,'linewidth',linewidth);   
        end
    end
end
