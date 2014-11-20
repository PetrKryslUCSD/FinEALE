function Commands =  targe2_mesher_vl_cmds(v)
    % Convert polyline  into targe2 commands.
    %
    % function Commands =  targe2_mesher_vl_cmds(v)
    %
    % Create  cell array of commands (strings) that define curves for triangulation of a
    % domain given as an array of vertices.
    %
    % Create a uniform triangulation of the domain given as an array of vertex
    % coordinates, v, one row per vertex. The vertices needs to be supplied in
    % counterclockwise order.
    %
    % Input: 
    % v= an array of vertex
    % coordinates, v, one row per vertex. The vertices needs to be supplied in
    % counterclockwise order.
    %
    % Output: 
    % Commands =  cell array of strings (commands)
    %
    nv =length(v);
    minSide=Inf;
    for i= 1:nv
        Commands{i}=sprintf('curve %d line %g %g %g %g\n',i,v(i,1),v(i,2),v(mod(i,nv)+1,1),v(mod(i,nv)+1,2));
        minSide=min([minSide,norm(v(i,:)-v(mod(i,nv)+1,:))]);
    end
end
