function [fens,fes] =T6_block(Length,Width,nL,nW,options)
% Mesh of a rectangle of T6 elements. 
% 
% function [fens,fes] =T6_block(Length,Width,nL,nW,options)
%  
% Examples: 
%     [fens,fes] = T6_block(3.1,2.2,3,2,[]);
%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
% If a different configuration of triangles is desired, run one of the
% mesh generation functions below and then convert  the mesh to T6
% triangles using T3_to_T6.
%
% See the visual gallery by running test_block
% 
% See also:  T3_blocku, T3_cblock,  T3_cblock2d, T3_crossblock, T3_ablock,
%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn, T3_to_T6

        [fens,fes] = T3_block(Length,Width,nL,nW,options);
        [fens,fes] = T3_to_T6(fens,fes,options);
    end