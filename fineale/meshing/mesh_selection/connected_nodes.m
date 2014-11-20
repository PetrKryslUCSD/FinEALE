% Extract the node numbers of the nodes  connected by given finite elements. 
%
% function cn = connected_nodes(fes)
%
% Extract the list of unique node numbers for the nodes that 
%  are connected by the finite element set fes. Note that it is assumed 
%  that all the FEs are of the same type (the same number of connected nodes
%  by each cell).
% 
% See also: connected_fes

function cn = connected_nodes(fes)
    conn= fes.conn;
    all_nodes = conn;
    cn = unique(reshape(all_nodes,prod(size(conn)), 1),'rows');
end
