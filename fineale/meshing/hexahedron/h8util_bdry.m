function f = h8util_bdry(t)
% Extract the boundary connectivity from a H8 connectivity.
%
% function f = h8util_bdry(t)
 %
 % t= connectivity of the hexahedra
% 
% Output:
% f = connectivity of the quadrilateral faces of the boundary, one per row
     % Form all hyperfaces, non-duplicates are boundary cells
    hypf=  [t(:,[1, 4, 3, 2]);t(:,[1, 2, 6, 5]);t(:,[2, 3, 7, 6]);t(:,[3, 4, 8, 7]);t(:,[4, 1, 5, 8]);t(:,[6, 7, 8, 5])];
    sorted_hypf = sort(hypf,2);
    %index vectors m and n such that b =sorted_hypf(m) and sorted_hypf = b(n)
    [b,m,n] = unique(sorted_hypf,'rows');
    u = find(histc(n,1:max(n))==1);
    f = hypf(m(u),:);
end
