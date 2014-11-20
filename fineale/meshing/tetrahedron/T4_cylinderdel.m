function [fens,fes] = T4_cylinderdel(Length,Radius,nL,nR)
% Tetrahedral (T4) Delaunay Mesh of a cylinder.
%
% function [fens,fes] = T4_cylinderdel(Length,Radius,nL,nR)
%
% Tetrahedral (T4) Delaunay Mesh of a cylinder of Length and Radius.
% Divided into elements: nL (lengthwise), nR (along the radius).
%
% Examples: 
% [fens,fes] = T4_cylinderdel(3.1,1.7,5,2);
% figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on

    dR=Radius/(nR+1/2);
    xnnodes= 3*(2*(0:1:nR)+1);
    nnodes=(nL+1)*sum(xnnodes);
    x=zeros(nnodes,1);
    y=zeros(nnodes,1);
    z=zeros(nnodes,1);
    xs=zeros(nnodes,3);
    k = 1;
    for j = 0:1:nL
        xjiggle = Length/nL/10*j*(j-nL)/nL/nL*rand;
        for n= 0:1:nR
            r =dR*(1/2+n);
            rjiggle=Radius/nR/10*rand*(nR-n)/nR;
            r= r+rjiggle;
            dA =2*pi/xnnodes(n+1);
            dAjiggle =dA/10*rand*(nR-n)/nR;
            for m=1:xnnodes(n+1)
                x(k)=j*Length/nL+xjiggle;
                y(k)=r*cos(dAjiggle+dA*(m-1));
                z(k)=r*sin(dAjiggle+dA*(m-1));
                xs(k,:)=[x(k),y(k),z(k)];
                k=k+1;
            end
        end
    end
    if (verLessThan('matlab', '7.6'))
        T = delaunay3(x,y,z);% 
    else
        T = DelaunayTri(x,y,z);% replaced DELAUNAY3 that will be removed in a future release
    end
    
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
