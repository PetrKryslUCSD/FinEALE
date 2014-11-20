function show_field_as_markers(gv,options)
% Visualize a field as markers. 
%
% function show_field_as_markers(gv,options)
%
% attributes: 
% x= geometry field 
% u= vector field to be visualized 
% nl= list of nodes, if empty then all nodes are used
%
% options may have further attributes which are interpreted by
% draw_marker().
%
    x= options.x;
    u= options.u;
    if isfield(options,'nl')
        nl = options.nl;
    else
        nl =(1: x.nfens);
    end
    p=gather_values(x,nl);
    v=gather_values(u,nl);
    for i= 1:length(nl)
        draw_marker (gv,p(i,:),options);
    end
end
