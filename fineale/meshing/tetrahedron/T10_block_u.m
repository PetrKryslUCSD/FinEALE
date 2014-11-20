function [fens,fes] = T10_block_u(Length,Width,Height,nL,nW,nH)
% Tetrahedral (T10) mesh of a rectangular block.
%
% function [fens,fes] = T10_block_u(Length,Width,Height,nL,nW,nH)
%
% Range =<0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z).
%
% See also: T4_blocka

    [fens,fes] = T4_blockdel(Length,Width,Height,nL,nW,nH);
    [fens,fes] = T4_to_T10(fens,fes);
end
