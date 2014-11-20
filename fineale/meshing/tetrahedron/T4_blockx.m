function [fens,fes] = T4_blockx(xs,ys,zs,orientation)
    % Graded tetrahedral (T4) mesh of a rectangular block.
    %
    % function [fens,fes] = T4_blockx(xs,ys,zs,orientation)
    %
    % xs,ys,zs =Locations of the individual planes of nodes.
    % orientation = either one of 'a', 'b', 'ca', 'cb' (when omitted, default is 'a')
    %
    % The mesh is produced by splitting each logical  rectangular cell into six
    % tetrahedra.
    %
    % Output:
    % fens= finite element node set
    % fes = finite element set
    %
    % Examples:
    %     [fens,fes] = T4_blockx(1/125*(0:1:5).^3,4+(0:2:8),[3, 5, 9]);
    %     drawmesh ({fens,fes})
    %
    %     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'a');
    %     figure; drawmesh({fens,fes},'fes','facecolor','red'); hold on
    %
    %     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'b');
    %     figure; drawmesh({fens,fes},'fes','facecolor','g'); hold on
    %
    %     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'ca');
    %     figure; drawmesh({fens,fes},'fes','facecolor','b'); hold on
    %
    %     [fens,fes] = T4_blockx([1,2,3,4].^2,[2,3,4,5].^2,[1,2,3,5].^2,'cb');
    %     figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
    %
    % See also: T4_blocka,  T4_blockb,  T4_blockca,  T4_blockcb
    if (~exist( 'orientation','var'))
        orientation ='a';
    end
    nL =length(xs)-1;
    nW =length(ys)-1;
    nH =length(zs)-1;
    nnodes=(nL+1)*(nW+1)*(nH+1);
    ncells=6*(nL)*(nW)*(nH);
    xyzs=zeros(nnodes,3);
    conns=zeros(ncells,4);
    if (strcmp(orientation,'a'))
        t4ia = [1, 8, 5, 6; 3, 4, 2, 7; 7, 2, 6, 8; 4, 7, 8, 2; 2, 1, 6, 8; 4, 8, 1, 2];
        t4ib = [1, 8, 5, 6; 3, 4, 2, 7; 7, 2, 6, 8; 4, 7, 8, 2; 2, 1, 6, 8; 4, 8, 1, 2];
    elseif (strcmp(orientation,'b'))
        t4ia = [2,7,5,6; 1,8,5,7; 1,3,4,8; 2,1,5,7; 1,2,3,7; 3,7,8,1];
        t4ib = [2,7,5,6; 1,8,5,7; 1,3,4,8; 2,1,5,7; 1,2,3,7; 3,7,8,1];
    elseif (strcmp(orientation,'ca'))
        t4ia = [8, 4, 7, 5; 6, 7, 2, 5; 3, 4, 2, 7; 1, 2, 4, 5; 7, 4, 2, 5];
        t4ib = [7, 3, 6, 8; 5, 8, 6, 1; 2, 3, 1, 6; 4, 1, 3, 8; 6, 3, 1, 8];
    elseif (strcmp(orientation,'cb'))
        t4ia = [7, 3, 6, 8; 5, 8, 6, 1; 2, 3, 1, 6; 4, 1, 3, 8; 6, 3, 1, 8];
        t4ib = [8, 4, 7, 5; 6, 7, 2, 5; 3, 4, 2, 7; 1, 2, 4, 5; 7, 4, 2, 5];
    end
    f=1;
    for k=1:(nH+1)
        for j=1:(nW+1)
            for i=1:(nL+1)
                xyzs(f,:)=[xs(i) ys(j) zs(k)];
                f=f+1;
            end
        end
    end

    fens=fenode_set(struct('xyz',xyzs));
    
    gc=1;
    for i=1:nL
        for j=1:nW
            for k=1:nH
                nn=node_numbers(i,j,k,nL,nW,nH);
                if (mod (sum( [i,j,k] ),2)==0)
                    t4i =t4ib;
                else
                    t4i =t4ia;
                end
                for r=1:size(t4i,1)
                    conns(gc,:)=nn(t4i(r,:));
                    gc=gc+1;
                end
            end
        end
    end
    fes =fe_set_T4(struct ('conn',conns(1:gc-1,:)));
end

function nn=node_numbers(i,j,k,nL,nW,nH)
    f=(k-1)*((nL+1)*(nW+1))+(j-1)*(nL+1)+i;
    nn=[f (f+1)  f+(nL+1)+1 f+(nL+1)];
    nn=[nn nn+((nL+1)*(nW+1))];
end
