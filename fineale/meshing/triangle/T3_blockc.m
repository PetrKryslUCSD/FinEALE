function [fens,fes] = T3_blockc(Length,Width,nL,nW,options)
% Mesh of a rectangle. Corners are bisected.
%
% function [fens,fes] = T3_blockc(Length,Width,nL,nW,options)
%
% Range =<0,Length> x <0,Width>
% Divided into triangular (t3) elements: nL, nW in the first, second (x,y).
% This function always bisects the corners.
% options = Attributes recognized by the constructor of fe_set_T3.
% 
% Examples: 
%     [fens,fes] = T3_blockc(3.1,2.2,3,2,[]);
%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
% See the visual gallery by running test_block
% 
% See also:  T3_blocku, T3_cblock,  T3_crossblock, T3_ablock,
%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn

    nnodes=(nL+1)*(nW+1);
    ncells=2*(nL)*(nW);
    xs=zeros(nnodes,2);
    conns=zeros(ncells,3);
    if ~isstruct(options)
        other_dimension = options; clear options;
        options.other_dimension = other_dimension;
    end
    f=1;
    for j=1:(nW+1)
        for i=1:(nL+1)
            xyz=[(i-1)*Length/nL (j-1)*Width/nW];
            xs(f,:)=xyz;
            f=f+1;
        end
    end
    fens=fenode_set(struct('xyz',xs));

   gc=1;
    for i=1:nL
        for j=1:nW
            if ((i==1) & (j==1)) | ((i==nL) & (j==nW))
                nn=node_numbers1(i,j,nL,nW);
                conns(gc,:)=nn;
                gc=gc+1;
                nn=node_numbers2(i,j,nL,nW);
                conns(gc,:)=nn;
                gc=gc+1;
            else
                nn=node_numbers3(i,j,nL,nW);
                conns(gc,:)=nn;
                gc=gc+1;
                nn=node_numbers4(i,j,nL,nW);
                conns(gc,:)=nn;
                gc=gc+1;
            end
        end
    end
    options.conn =conns;
    fes =fe_set_T3(options);
end

   
function nn=node_numbers1(i,j,nL,nW)
    f=(j-1)*(nL+1)+i;
    nn=[f (f+1) f+(nL+1)+1];
    return;
end
function nn=node_numbers2(i,j,nL,nW)
    f=(j-1)*(nL+1)+i;
    nn=[f  f+(nL+1)+1 f+(nL+1)];
    return;
end

function nn=node_numbers3(i,j,nL,nW)
    f=(j-1)*(nL+1)+i;
    nn=[f (f+1) f+(nL+1)];
    return;
end
function nn=node_numbers4(i,j,nL,nW)
    f=(j-1)*(nL+1)+i;
    nn=[(f+1)  f+(nL+1)+1 f+(nL+1)];
    return;
end
