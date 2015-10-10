function [fens,fes] = T10_blockcb(Length,Width,Height,nL,nW,nH)
% Tetrahedral (T10) mesh of a rectangular block.
%
% function [fens,fes] = T10_blockcb(Length,Width,Height,nL,nW,nH)
%
% Range =<0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z).
%
% See also: T4_blockcb

    [fens,fes] = T10_blockr(Length,Width,Height,nL,nW,nH,@T4_blockcb);
end
