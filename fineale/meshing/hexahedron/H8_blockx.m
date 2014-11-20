function [fens,fes] = H8_blockx(xs,ys,zs)
% Graded mesh of a 3-D block of H8 finite elements.
  %
% Mesh of a 3-D block, H8 cells. The nodes are located at the Cartesian
% product of the three intervals on the input.  This allows for
% construction of graded meshes.
%
% function [fens,fes] = H8_blockx(xs,ys,zs)
%
% xs,ys,zs =Locations of the individual planes of nodes. 
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% Examples:
% 
%  [fens,fes] = H8_blockx(1/125*(0:1:5).^3,4+(0:2:8),[3, 5, 9]); drawmesh ({fens,fes})
%  
%  [fens,fes] = H8_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2)
%  drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_block, H8_composite_plate

    nL =length(xs)-1;
    nW =length(ys)-1;
    nH =length(zs)-1;

    nnodes=(nL+1)*(nW+1)*(nH+1);
    xyz =zeros(nnodes,3);
    ncells=(nL)*(nW)*(nH);
    conns =zeros(ncells,8);

    f=1;
    for k=1:(nH+1)
        for j=1:(nW+1)
            for i=1:(nL+1)
                xyz(f,:)=[xs(i) ys(j) zs(k)];
                f=f+1;
            end
        end
    end

    gc=1;
    for i=1:nL
        for j=1:nW
            for k=1:nH
                nn=node_numbers(i,j,k,nL,nW,nH);
                conns(gc,:)=nn;
                gc=gc+1;
            end
        end
    end
    fens=fenode_set(struct('xyz',xyz));
    fes =fe_set_H8(struct ('conn',conns));
end

% first go along Length, then Width, and finally Height
function nn=node_numbers(i,j,k,nL,nW,nH)
    f=(k-1)*((nL+1)*(nW+1))+(j-1)*(nL+1)+i;
    nn=[f (f+1)  f+(nL+1)+1 f+(nL+1)];
    nn=[nn nn+((nL+1)*(nW+1))];
end
