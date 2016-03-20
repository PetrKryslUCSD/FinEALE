function draw_colorbar(self, context)
% Draw the color bar
%
% function draw_colorbar(self, context)
%
% context= Structure with optional fields:
%    colormap=Colormap
%    position=where to put it in the normalized coordinates
%    minmax=minimum and maximum (also can be supplied as a field range)
%    label= descriptive text
%    fontname= name of the font to use for labels
%    fontsize= size of the fonts to use the labels
%    interpreter=interpreter to use for the label text (default is to let the
%      handle graphics decide; alternative is 'latex') 
    cmap=colormap;
    if isfield(context,'colormap')
        cmap = context.colormap;
    end
    cmap =cmap(2:end-1,:);
    color=[0, 0, 0];
    if isfield(context,'color')
        color = context.color;
    end
    FontName= get(0,'DefaultAxesFontName' );
    if isfield(context,'fontname')
        FontName = context.fontname;
    end
    Fontsize=12;
    if isfield(context,'fontsize')
        Fontsize = context.fontsize;
    end
    Interpreter=[];
    if isfield(context,'interpreter')
        Interpreter = context.interpreter;
    end
    position =[0.869 0.15 0.05 0.7];
    if isfield(context,'position')
        position = context.position;
    end
    minmax = [0, 1];
    if isfield(context,'minmax')
        minmax= context.minmax;
    end
    if isfield(context,'range')
        minmax= context.range;
    end
    if (min(minmax)==max(minmax))% If the range is empty, enlarge it a bit down and up
        minmax=[min(minmax)-eps(min(minmax)),max(minmax)+eps(max(minmax))];
    end
    label = ' Data';
    if isfield(context,'label')
        label= context.label;
    end
    colormap(cmap);
    cbh=colorbar;
    set(cbh,...
        'Position',position,'YLim',minmax,'YTickMode','auto');
    YTick= get(cbh,'YTick');
    YTick=(YTick-minmax(1))./(minmax(2)-minmax(1));
    YTickLabel= get(cbh,'YTickLabel');
    m=size(YTickLabel,1);
    yl(1:m)=deal({' '});
    for i=1:m
        yl(i)={num2str(YTick(i)*(minmax(2)-minmax(1))+minmax(1))};
    end
    YTickLabel=yl;
    mtd=1/5/length(YTick);
    if abs(YTick(1)) >mtd
        YTick=[0,YTick];
        YTickLabel=cat(2,{num2str(minmax(1))},YTickLabel);
    end
    if abs(YTick(end)-1) >mtd
        YTick=[YTick,1];
        YTickLabel= cat(2,YTickLabel,{num2str(minmax(2))});
    end
    XTickLabel= get(cbh,'XTickLabel');
    m=size(XTickLabel,1);
    xl(1:m)=deal({' '});
    %     set(cbh,...
    %         'Position',position,'XColor',color,'YColor',color,...
    %         'YLim',[0, 1],'YTick',YTick,'YTickLabel',YTickLabel,...
    %         'XTickMode','manual','XTickLabelMode','manual','XTick',[],'XTickLabel',xl);
    set(cbh,...
        'Position',position,'XColor',color,'YColor',color,...
        'Limits',[0, 1]);
    set(cbh,...
        'Position',position,'XColor',color,'YColor',color,...
        'Ticks',YTick,'TickLabels',YTickLabel);
    try% R2014b deleted this property
        set(cbh,'CLim',[0, 1]);
    catch
    end
    if (~isempty(Fontsize))
        set(cbh,'Fontsize',Fontsize);
    end
    if (~isempty(FontName))
        set(cbh,'fontname',FontName);
    end
    set(get(cbh,'XLabel'),'String',label);
    if (~isempty(Fontsize))
        set(get(cbh,'XLabel'),'Fontsize',Fontsize);
    end
    if (~isempty(FontName))
        set(get(cbh,'XLabel'),'FontName',FontName);
    end
    if (~isempty(Interpreter))
        set(get(cbh,'XLabel'),'Interpreter',Interpreter);
    end
end
