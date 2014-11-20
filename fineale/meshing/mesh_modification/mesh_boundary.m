function bdry_fes = mesh_boundary (fes, options)
% Extract the boundary finite elements from a mesh.
%
% function bdry_fes = mesh_boundary(fes, options)
%
% Extract the finite elements of manifold dimension (n-1) from the
% supplied list of finite elements of manifold dimension (n).
%    options = struct with any attributes that should be passed to the
%    constructor of the boundary finite elements
% 
%
    if (~exist('options','var'))
        options =[];
    end
    make =fes.boundary_fe();
    % Form all hyperfaces, non-duplicates are boundary cells
    hypf= fes.boundary_conn();
    sorted_hypf = sort(hypf,2);
    %index vectors m and n such that b =sorted_hypf(m) and sorted_hypf = b(n)
    [b,m,n] = unique(sorted_hypf,'rows');
    u = find(histc(n,1:max(n))==1);
    bdry = hypf(m(u),:);
    nbc= size(bdry, 1);
    if (~isempty(options))
        Parameters = options;
    end
    Parameters.conn =bdry;
    bdry_fes=make(Parameters);
    return
end
