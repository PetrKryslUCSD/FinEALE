function draw_annotation(self, xywh, txt, context)
% Draw annotation text.
%
% function draw_annotation(self, xy, txt, context)
%
%
% context = struct
% with mandatory fields
% xy = anchor location and width and height (in normalized coordinates)
% txt = text
% context = struct
% with optional fields
% color = color of the text
% offset = offset with respect to the text anchor location xyz
%
    color ='black';
    if isfield(context,'color')
        color=context.color;
    end
    offset = 0;
    if isfield(context,'offset')
        offset=context.offset;
    end
    fontName = [];
    if isfield(context,'fontname')
        fontName=context.fontname;
    end
    h=annotation('textbox',xywh, 'String',txt,'Color', color);
    if isfield(context,'backgroundcolor')
       set(h,'backgroundcolor',context.backgroundcolor);
    end
    if (~isempty(fontName))
        set(h,'fontname',fontName);
    end
    if isfield(context,'fontsize')
       set(h,'fontsize',context.fontsize);
    end
    if isfield(context,'fontangle')
       set(h,'fontangle',context.fontangle);
    end
    if isfield(context,'interpreter')
       set(h,'interpreter',context.interpreter);
    end
    return;
end