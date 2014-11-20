function draw_box(self, box, context)
% Draw Box.
%
% function draw_box(self, box, context)
%
% required arguments:
% box=[min(x(:,1)),max(x(:,1))], or
% box=[min(x(:,1)),max(x(:,1)),min(x(:,2)),max(x(:,2))], or
% box=[min(x(:,1)),max(x(:,1)),min(x(:,2)),max(x(:,2)),min(x(:,3)),max(x(:,3))] 
% context = struct with following optional fields
%    edgecolor = color if polyline should be drawn in solid color,
%    linewidth = line width, integer
%    length_units= number, representing the length units in which the
%         graphic should be displayed
%
    edgecolor='black';
    if isfield(context,'edgecolor')
        edgecolor =context.edgecolor;
    end
    linewidth=1;
    if isfield(context,'linewidth')
        linewidth =context.linewidth;
    end
    if isfield(context,'length_units')
        box=box/context.length_units;
    end
    if (~strcmp(edgecolor,'none'))
        if (length(box)==6)
            X =  box([1,2,2,1,1]);
            Y =  box([3,3,4,4,3]);
            Z =  box([5,5,5,5,5]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            Z =  box([6,6,6,6,6]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            X =  box([1,1]);
            Y =  box([3,3]);
            Z =  box([5,6]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            X =  box([2,2]);
            Y =  box([3,3]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            Y =  box([4,4]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
            X =  box([1,1]);
            line(X',Y',Z', 'Color',edgecolor,'linewidth',linewidth);
        elseif (length(box)==4)
            X =  reshape(x(edge,1),size(edge,1),size(edge,2));
            Y =  reshape(x(edge,2),size(edge,1),size(edge,2));
            line(X',Y', 'Color',edgecolor,'linewidth',linewidth);
        end
    end
end
