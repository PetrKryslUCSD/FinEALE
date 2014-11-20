function [fens,fesa] = merge_n_meshes(fensa, fesa, tolerance)
    % Merge several meshes together.
%
% function [fens,fesa] = merge_n_meshes(fensa, fesa, tolerance)
%
% Merge several meshes together either by simple concatenation of nodes or by 
% gluing together nodes within tolerance. 
%
% Inputs:
% fensa= cell array of node sets, one for each mesh;
% fesa= cell array of finite element sets, one for each mesh;  
% tolerance= Geometric tolerance, maybe supplied as zero (>=0). 
% 
% The meshes are glued together by
% merging the nodes that fall within a box of size "tolerance". If tolerance 
% is set to zero, no merging of nodes is performed; the nodes from the meshes are 
% simply concatenated together.
% 
% The merged node set, fens, and the cell array of finite element sets with
% renumbered  connectivities are returned.
% 
% Outputs:
% fens= merged node set, 
% fesa= cell array of finite element sets updated to use the merged node set.
%
%
% See also: merge_meshes


    if (length(fensa))~=(length(fesa))
        error('(length(fensa))~=(length(fesa))');
    end
    if (length(fensa))==1
        fens=fensa{1};
        return
    end
    fens=fensa{1};
    for j=2:length(fesa)
        [fens,new_indexes_of_fens1_nodes] = fuse_nodes(fensa{j}, fens, tolerance);
        fesa{j}= update_conn(fesa{j},new_indexes_of_fens1_nodes);
    end
end
