function [fens,fes] = T4_delx(xyz)
% Tetrahedral (T4) Delaunay Mesh of the supplied locations of vertices.
%
% function [fens,fes] = T4_delx(xyz)
%
% Tetrahedral (T4) Delaunay Mesh of a set of vertices.
% Divided into elements: xyz (locations of vertices, one per row), free
% (array that indicates in which directions the nodes are free to move, 1,
% and in which they are not allowed to move,  0).  
%
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
% Unstructured mesh of a rectangular block.
%       nL=4; nW=3; nH=7;
%       [fens,~] = H8_block(nL,nW,nH,nL,nW,nH);
%       free= ones(size(fens.xyz));
%       nl=fenode_select(fens,struct( 'box', [0,0,-inf,+inf,-inf,+inf],'inflate',1/nL/10));
%       free(nl,1)=0;
%       nl=fenode_select(fens,struct( 'box', [nL,nL,-inf,+inf,-inf,+inf],'inflate',1/nL/10));
%       free(nl,1)=0;
%       nl=fenode_select(fens,struct( 'box', [-inf,+inf,0,0,-inf,+inf],'inflate',1/nW/10));
%       free(nl,2)=0;
%       nl=fenode_select(fens,struct( 'box', [-inf,+inf,nW,nW,-inf,+inf],'inflate',1/nW/10));
%       free(nl,2)=0;
%       nl=fenode_select(fens,struct( 'box', [-inf,+inf,-inf,+inf,0,0],'inflate',1/nH/10));
%       free(nl,3)=0;
%       nl=fenode_select(fens,struct( 'box', [-inf,+inf,-inf,+inf,nH,nH],'inflate',1/nH/10));
%       free(nl,3)=0;
%       r=(rand(size(free))- 0.5)*2;
%       fens.xyz=fens.xyz+r.*free*0.1;
%       [fens,fes] = T4_delx(fens.xyz);
%       drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: T4_blockdel
 
  box = bounding_box(xyz);
  tol=min([box(2)-box(1),box(4)-box(3),box(6)-box(5)])/size(xyz,1)/10;

fens=fenode_set(struct('xyz', xyz));
xs=fens.xyz;
T = DelaunayTri(xs(:,1),xs(:,2),xs(:,3));% replaced DELAUNAY3 that will be removed in a future release

ncells=size(T,1);
conns=zeros(ncells,4);
for i=1:ncells
    if tetvol(xs(T(i,:),:)) > 0
        conns(i,:) =T(i,:);
    else
        conns(i,:) =T(i,[1, 3, 2, 4]);
    end
end
fens.xyz=xs;
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
