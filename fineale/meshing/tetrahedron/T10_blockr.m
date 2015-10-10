function [fens,fes] = T10_blockr(Length,Width,Height,nL,nW,nH,hT4_block)
% Tetrahedral (T10) mesh of a rectangular block.
%
% function [fens,fes] = T10_blockr(Length,Width,Height,nL,nW,nH)
%
% Range =<0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z).
%
% See also: T10_blocka, T10_blockb, T10_blockca, T10_blockcb

    [fens,fes] = hT4_block(Length,Width,Height,nL,nW,nH);
    [fens,fes] = T4_to_T10(fens,fes);
end
