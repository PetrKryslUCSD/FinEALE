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

% Approach 1
% sorted_hypf = sort(hypf,2);
% %index vectors m and n such that b =sorted_hypf(m) and sorted_hypf = b(n)
% [b,m,n] = unique(sorted_hypf,'rows');
% u = find(histc(n,1:max(n))==1);
% bdryconn = hypf(m(u),:);

% Approach 2
bdryconn =myunique(hypf);

if (~isempty(options))
    Parameters = options;
end
Parameters.conn =bdryconn;
bdry_fes=make(Parameters);
return
end

function  Out  =myunique(A)
    sA=sort(A,2);
    [sA,rix]  = sortrows(sA);;
    d=sA(1:end-1,:)~=sA(2:end,:);
    ad=[true; any(d,2)];
    iu =find((ad&[ad(2:end);true])==true);
    Out =A(rix(iu),:);
end
