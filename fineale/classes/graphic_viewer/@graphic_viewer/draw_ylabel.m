function draw_ylabel(self, lbl, context)
% Draw Y label.
%
% function draw_ylabel(self, lbl, context)
%
% required arguments:
% lbl=label 
% context = struct with following optional fields
%    color = color,
%    linewidth = line width, integer
%
    fontName=[];
    if (self.pk_defaults )
        fontName='Times';
    end
    l=ylabel(lbl);
    if (~isempty(fontName))
        set(l,'fontname',fontName);
    end
    if isfield(context,'fontsize')
       set(l,'fontsize',context.fontsize);
    end
    if isfield(context,'fontangle')
       set(l,'fontangle',context.fontangle);
    end
    if isfield(context,'interpreter')
        set(l,'interpreter',context.interpreter);
    else
        set(l,'interpreter','latex');
    end
end
