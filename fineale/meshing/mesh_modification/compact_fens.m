function [fens, new_numbering] =compact_fens(fens, connected)
% Compact the finite element node set.
%
% Compact the finite element node set by deleting unconnected nodes.
% 
% function [fens, new_numbering] =compact_fens(fens, connected)
% 
% fens = array of finite element nodes
% connected = The array element connected(j) is either 0 (when j is an unconnected
%    node), or a positive number (when node j is connected to other nodes by
%    at least one finite element)  
%
% Output:
% fens = new array of finite element nodes
% new_numbering= array which tells where in the new fens array the
% connected nodes are (or 0 when the node was unconnected). For instance,
% node 5 was connected, and in the new array it is the third node: then
% new_numbering(5) is 3.
% 
% Examples: 
%
% Let us say there are nodes not connected to any finite element that you
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
    nconnfens =find(connected~=0);
    new_numbering=zeros(count(fens),1);
    xyz =fens.xyz;
    nid =(1:count(fens))';
    nxyz =xyz;
    id=1;
    for i=1:length(connected)
        if(connected(i)),
            new_numbering(i)=id; 
            nxyz(id,:) =xyz (i,:);
            nid(id) =id;
            id=id+1;
        end,
    end
    fens=fenode_set (struct('xyz',nxyz(1:id-1,:)));
end
