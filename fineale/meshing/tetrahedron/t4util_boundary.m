function f = t4util_boundary(t)
% Extract the boundary connectivity from a T4 connectivity.
%
% function f = t4util_boundary(t)
%
% t= connectivity of the tetrahedra
%
% Output:
% f  = connectivity of the boundary faces
%
% See also:

     % Form all hyperfaces, non-duplicates are boundary cells
    hypf=  [t(:,[1, 3, 2]);t(:,[1, 2, 4]);t(:,[2, 3, 4]);t(:,[1, 4, 3])];
    sorted_hypf = sort(hypf,2);
    %index vectors m and n such that b =sorted_hypf(m) and sorted_hypf = b(n)
    [b,m,n] = unique(sorted_hypf,'rows');
    u = find(histc(n,1:max(n))==1);
    f = hypf(m(u),:);
end
