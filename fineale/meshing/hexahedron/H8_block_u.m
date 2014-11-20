function [fens,fes] = H8_block_u(Length,Width,Height,nL,nW,nH)
% Mesh of a block, unstructured, obtained by splitting of tetrahedra.
%
% function [fens,fes] = H8_block_u(Length,Width,Height,nL,nW,nH)
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
%     [fens,fes] = H8_block_u(2,3,4, 2,6,4);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_block, H8_blockx

    nL=max([round(nL/2),1]);
    nW=max([round(nW/2),1]);
    nH=max([round(nH/2),1]);
    [fens,fes] = T4_blockdel(Length,Width,Height,nL,nW,nH);
    [fens,fes] = T4_to_H8(fens,fes);
end
