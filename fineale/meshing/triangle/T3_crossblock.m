function [fens,fes] = T3_crossblock(Length,Width,nL,nW,options)
% Mesh of a rectangle located in the range (cross-diagonal).
%
% function [fens,fes] = T3_crossblock(Length,Width,nL,nW,options)
%
% <0,Length> x <0,Width>
% Divided into triangular (t3) elements: nL, nW in the first, second (x,y).
% options = Attributes recognized by the constructor of fe_T3.
% 
% Examples: 
%     [fens,fes] = T3_crossblock(3.1,2.2,3,2,[]);
%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
% See the visual gallery by running test_block
% 
% See also:  T3_blocku, T3_cblock,   T3_crossblock, T3_ablock,
%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn

    if ~isstruct(options)
        other_dimension = options; clear options;
        options.other_dimension = other_dimension;
    end
    tol=0.01/max( [nL,nW]);
    fens= fenode_set(struct('xyz',[[0 0];[1 0];[1 1];[0 1];[1/2 1/2]]));
    options.conn =[[1, 2, 5];[2, 3, 5];[3, 4, 5];[4, 1, 5]];
    fes=fe_set_T3(options);
    fens1 =fens;
    fes1=fes;
    for i= 1:nL-1
        fes2=fes1;
        fens2 = transform_apply(fens1,@translate1, i);
        [fens,fes,fes2] = merge_meshes(fens, fes, fens2, fes2, tol);
        fes= cat(fes,fes2);
    end
    fens1 =fens;
    fes1=fes;
    for i= 1:nW-1
        fes2=fes1;
        fens2 = transform_apply(fens1,@translate2, i);
        [fens,fes,fes2] = merge_meshes(fens, fes, fens2, fes2, tol);
        fes= cat(fes,fes2);
    end
    fens = transform_apply(fens,@scale, [Length/nL,Width/nW]);
end

function xyz=translate1(xyz, i)
    xyz= [xyz(1)+i,xyz(2)];
end

function xyz=translate2(xyz, i)
    xyz= [xyz(1),xyz(2)+i];
end


function xyz=scale(xyz, scales)
    xyz= [xyz(1)*scales(1),xyz(2)*scales(2)];
end

