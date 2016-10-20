function [fens,fes] = T4_blockdelx(Length,Width,Height,nL,nW,nH,shiftamplitude)
% Tetrahedral (T4) Delaunay Mesh of a rectangular block.
%
% function [fens,fes] = T4_blockdelx(Length,Width,Height,nL,nW,nH,shiftamplitude)
%
% Tetrahedral (T4) Delaunay Mesh of a cylinder of Length and Radius.
% Divided into elements: nL (lengthwise), nR (along the radius).
%
% Mesh of a block located in a given range, H8 cells
% <0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z). H8 elements obtained by subdividing tetrahedra.
% shiftamplitude  =  amplitude of the random node  shift;  default is 1/10
%      of the node spacing
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples:
% With the default node perturbations:
%     [fens,fes] = T4_blockdelx(2,3,4, 3,6,4);
%       drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% With specified amplitude of the perturbation of the node locations: 
%    [fens,fes] = T4_blockdelx(2,3,4, 3,6,4, 0.7);
%       drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: T4_blockdel, perturb_nodes

% First generate  a dummy mesh  just to get the locations of the nodes
[fens,fes] = H8_block(nL,nW,nH,nL,nW,nH);
tol=1/100;

if (~exist('shiftamplitude', 'var'))
    shiftamplitude=10*tol;
end

free= ones(size(fens.xyz));
xs=fens.xyz;
for i=1:size(xs,1)
    if ~((abs(xs(i,1))>tol)&&(abs(xs(i,1)-nL)>tol))
        %             Node not movable  in the X direction
        free(i,1)=0;
    end
    if ~((abs(xs(i,2))>tol)&&(abs(xs(i,2)-nW)>tol))
        %             Node not movable  in the Y direction
        free(i,2)=0;
    end
    if ~((abs(xs(i,3))>tol)&&(abs(xs(i,3)-nH)>tol))
        %             Node not movable  in the Z direction
        free(i,3)=0;
    end
end

[fens] = perturb_nodes(fens, shiftamplitude*[1,1,1], free);
xs=fens.xyz;
T = DelaunayTri(xs(:,1),xs(:,2),xs(:,3));% replaced DELAUNAY3 that will be removed in a future release

xs=[Length/nL*xs(:,1),Width/nW*xs(:,2),Height/nH*xs(:,3)];

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
