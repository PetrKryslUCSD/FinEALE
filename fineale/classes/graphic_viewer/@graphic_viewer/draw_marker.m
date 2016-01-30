function draw_marker(self, c, context)
  % Draw marker.
  %
  % function draw_marker(self, c,  context)
  %
  % required arguments:
  % c= coordinates of the center,  one point per row;
  % context = struct with following optional fields
  %    markersize= marker size;
  %    color = color if marker should be drawn in solid color
  %    length_units= number, representing the length units in which the
  %         graphic should be displayed
  %
  dim=size(c,2);
  if dim<3
    c(:,end+1)=0;
  end
  if dim<2
    c(:,end+1)=0;
  end
  if isfield(context,'length_units')
    c=c/context.length_units;
  end
  if isfield(context,'color')
    color= context.color;
  else
    color='k';
  end
  if isfield(context,'marker')
    marker= context.marker;
  else
    marker='o';
  end
  if isfield(context,'markersize')
    markersize= context.markersize;
  else
    markersize=3;
  end
  if isfield(context,'linewidth')
    linewidth= context.linewidth;
  else
    linewidth=1;
  end
  line(c(:,1),c(:,2),c(:,3),'linestyle','none','linewidth',linewidth,'color',color,'marker',marker,'markersize',markersize);
end

