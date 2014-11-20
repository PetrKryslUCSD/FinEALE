function [fens,fes] = T4_blockdel(Length,Width,Height,nL,nW,nH)
% Tetrahedral (T4) Delaunay Mesh of a rectangular block.
%
% function [fens,fes] = T4_blockdel(Length,Width,Height,nL,nW,nH)
%
% Tetrahedral (T4) Delaunay Mesh of a cylinder of Length and Radius.
% Divided into elements: nL (lengthwise), nR (along the radius).
%
% Mesh of a block located in a given range, H8 cells
% <0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z). H8 elements obtained by subdividing tetrahedra.
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes] = T4_blockdel(2,3,4, 2,6,4);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also:

tol=1/100;
mtol=80*tol;

    [fens,fes] = H8_block(nL,nW,nH,nL,nW,nH);
    xs=fens.xyz;
    for i=1:size(xs,1)
        if (abs(xs(i,1))>tol)&&(abs(xs(i,1)-nL)>tol)
            xs(i,1)=xs(i,1)+mtol*(rand(1,1)-0.5);
        end
        if (abs(xs(i,2))>tol)&&(abs(xs(i,2)-nW)>tol)
            xs(i,2)=xs(i,2)+mtol*(rand(1,1)-0.5);
        end
        if (abs(xs(i,3))>tol)&&(abs(xs(i,3)-nH)>tol)
            xs(i,3)=xs(i,3)+mtol*(rand(1,1)-0.5);
        end
            
    end
    fens.xyz=xs;
    
    if (verLessThan('matlab', '7.6'))
        T = delaunay3(xs(:,1),xs(:,2),xs(:,3));% 
    else
        T = DelaunayTri(xs(:,1),xs(:,2),xs(:,3));% replaced DELAUNAY3 that will be removed in a future release
    end
    
    xs=[Length/nL*xs(:,1),Width/nW*xs(:,2),Height/nH*xs(:,3)];
    
    ncells=size(T,1);
    conns=zeros(ncells,4);
    for i=1:ncells
        xyz=zeros (4, 3);
        for k=1:4
            xyz(k,:) = xs(T(i,k),:);
        end
        if tetvol(xyz) > 0
            conns(i,:) =T(i,:);
        else
            conns(i,:) =T(i,[1, 3, 2, 4]);
        end
%         tetvol(xyz)
    end
    fens=fenode_set(struct('xyz',xs));
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
