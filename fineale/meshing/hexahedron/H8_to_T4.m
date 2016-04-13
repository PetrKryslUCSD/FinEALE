function [fens,fes] = H8_to_T4(fens,fes,sxy,sxz,syz)
% Convert a mesh of hexahedra H8 to tetrahedra T4.
%
% function [fens,fes] = H8_to_T4(fens,fes)
%
% Arguments and
% Output:
% fens= finite element node set
% fes = finite element set
% sxy,sxz,syz=amount of shear to be applied to make  it possible to
%   triangulate a regular arrangement of nodes with the Delaunay
%   triangulation.
% 
% Note: the labels of the original hexahedra are not transferred to the resulting tetrahedra.
M=[0,sxy,sxz; sxy,0,syz; sxz,syz,0];
X=fens.xyz;
X(:,:)=X(:,:)+X(:,:)*M;
if (verLessThan('matlab', '7.6'))
    T = delaunay3(X(:,1),X(:,2),X(:,3));%
else
    T = DelaunayTri(X(:,1),X(:,2),X(:,3));% replaced DELAUNAY3 that will be removed in a future release
end
% X=fens.xyz;
ncells=size(T,1);
conns=zeros(ncells,4);
for i=1:ncells
    xyz=zeros (4, 3);
    for k=1:4
        xyz(k,:) = X(T(i,k),:);
    end
    if tetvol(xyz) > 0
        conns(i,:) =T(i,:);
    else
        conns(i,:) =T(i,[1, 3, 2, 4]);
    end
end
fens=fenode_set(struct('xyz',X));
fes =fe_set_T4(struct ('conn',conns));
end

%compute volume of a tetrahedron
% Given the 4x3 vertex coordinate matrix V of a tetrahedron, TETVOL(V)
% returns the volume of the tetrahedron.
function vol = tetvol(v)
    vol = det([v(2,:)-v(1,:);v(3,:)-v(1,:);v(4,:)-v(1,:)])/6;
    %     if abs (vol) < 0.1
    %         warning (' sliver?')
    %     end
    return;
end
