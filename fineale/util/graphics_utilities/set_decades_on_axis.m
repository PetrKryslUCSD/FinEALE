function  set_decades_on_axis (ax,which)
if (~exist('which','var'))
    which ='all';
end

if (strcmp(which,'all' ))||(strcmp(which,'x' ))
    xlim=get(ax,'xlim');
    xlim=10.^[floor(log10(xlim(1))),ceil(log10(xlim(2)))];
    set(ax,'xlim',xlim);
    set(ax,'xtick',10.^[floor(log10(xlim(1))):ceil(log10(xlim(2)))]);
end

if (strcmp(which,'all' ))||(strcmp(which,'y' ))
    ylim=get(ax,'ylim');
    ylim=10.^[floor(log10(ylim(1))),ceil(log10(ylim(2)))];
    set(ax,'ylim',ylim);
    set(ax,'ytick',10.^[floor(log10(ylim(1))):ceil(log10(ylim(2)))]);
end

