% Renumber the nodes in the connectivity of the finite elements based on a new
% numbering for the nodes.
% 
% function fes = renumber_fe_conn(fes, new_numbering)
% 
% fes =finite element set 
% new_numbering = new serial numbers for the nodes.  The connectivity
%           should be changed as conn(j) --> new_numbering(conn(j))
% 
% Output:
% fes = new Finite element set
% 
% Let us say there are nodes not connected to any Finite element that you
% would like to remove from the mesh: here is how that would be
% accomplished.
% 
% connected = find_unconn_fens(fens, fes);
% [fens, new_numbering] =compact_fens(fens, connected);
% fes = renumber_fe_conn(fes, new_numbering);
% 
% Finally, check that the mesh is valid:
% validate_mesh(fens, fes);
% 
function fes = renumber_fe_conn(fes, new_numbering)
    conn =fes.conn;
    for i=1:size(conn,1)
        c=conn(i,:);
        conn(i,:) =new_numbering(c)';
    end
    fes.conn=conn;
end