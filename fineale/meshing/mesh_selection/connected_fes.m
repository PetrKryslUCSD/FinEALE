function cg = connected_fes(fes, node_list)
% Extract the numbers of the finite elements connected to given nodes. 
%
% function cg = connected_fes(fes, node_list)
%
% Extract the list of numbers for the fes  that 
%  are connected to given nodes. 
% 
% Warning: this tends to be an expensive operation.
% 
% See also: connected_nodes
%
conn= fes.conn;
cg=zeros(size(conn,1),1);
for j=1:size(conn,1)
	cg(j)= (~isempty( intersect(conn(j,:), node_list) ));
end
cg =find(cg~=0);
end
