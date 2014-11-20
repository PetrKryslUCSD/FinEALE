function draw_text(self, xyz, txt, context)
% Draw text.
%
% function draw_text(self, xyz, txt, context)
%
%
% context = struct
% with mandatory fields
% xyz = anchor location
% txt = text
% context = struct
% with optional fields
%    color = color of the text
%    offset = offset with respect to the text anchor location xyz
%    length_units= number, representing the length units in which the
%         graphic should be displayed
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
    fontsize = [];
    if isfield(context,'fontsize')
        fontsize=context.fontsize;
    end
    annotate=0;
    if isfield(context,'annotate')
        annotate=context.annotate;
    end
    xyz=xyz + offset;
    xyz=[xyz 0 0];% just in case the node is in one- or two a dimensions
    if isfield(context,'length_units')
        xyz=xyz/context.length_units;
    end
    h=text(xyz(1),xyz(2),xyz(3),txt,'Color', color);
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