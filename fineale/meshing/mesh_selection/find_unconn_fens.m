function connected = find_unconn_fens(fens, fes)
% Find nodes that are not connected to any finite element.
% 
% function connected = find_unconn_fens(fens, fes)
% 
% fens = array of finite element nodes
% 
% Output:
% connected = array is returned which is for the node k either true (node k is
%      connected), or false (node k is not connected).
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
    connected=zeros(count(fens),1);
    fen2fe_map=fenode_to_fe_map (struct ('fes',fes));
    gmap=fen2fe_map.map;
    for i=1:length(gmap),
        connected(i) =(~isempty(gmap{i}));
    end
end
