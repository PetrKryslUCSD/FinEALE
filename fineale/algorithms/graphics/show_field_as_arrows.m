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
% draw_arrow().
%
    x= options.x;
    u= options.u;
    if isfield(options,'nl')
        nl = options.nl;
    else
        nl =(1: get (x,'nfens'));
    end
    p=gather(x,nl,'values','noreshape');
    v=gather(u,nl,'values','noreshape');
    for i= nl
        draw_arrow (gv,p(i,:),v(i,:),options);
    end
end
