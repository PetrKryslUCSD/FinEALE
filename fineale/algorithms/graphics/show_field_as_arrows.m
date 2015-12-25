function show_field_as_arrows(gv,options)
% Visualize a field as arrows. 
%
% function show_field_as_arrows(gv,options)
%
% attributes: 
% x= geometry field 
% u= vector field to be visualized 
% nl= list of nodes, if empty then all nodes are used
%
% options may have further attributes which are interpreted by
% the draw_arrow() method of the graphic viewer.
%
    x= options.x;
    u= options.u;
    if isfield(options,'nl')
        nl = options.nl;
    else
        nl =(1:x.nfens);
    end
    p=gather_values(x,nl);
    v=gather_values(u,nl);
    for i= nl
        draw_arrow (gv,p(i,:),v(i,:),options);
    end
end
