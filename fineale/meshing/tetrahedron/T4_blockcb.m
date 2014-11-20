function [fens,fes] = T4_blockcb(Length,Width,Height,nL,nW,nH)
    % Tetrahedral mesh of a rectangular block; orientation 'cb'.
    %
    % function [fens,fes] = T4_blockcb(Length,Width,Height,nL,nW,nH)
    %
    % The orientation of the tetrahedra in the logical rectangular sub
    % blocks is varied to produce a kind of cross-hatched pattern. 
    %
    % Range   =<0,Length> x <0,Width> x <0,Height> Divided into elements:
    % nL, nW, nH in the first, second, and third direction (x,y,z).
    %
    % Output:
    % fens= finite element node set
    % fes = finite element set
    %
    % Examples:
    %         [fens,fes] = T4_blockcb(2,3,4, 2,6,4);
    %         figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
    %
    % See also: T4_blockx,  T4_blocka,  T4_blockb,  T4_blockca
    [fens,fes] = T4_blockx(linspace(0,Length,nL+1)',linspace(0,Width,nW+1)',linspace(0,Height,nH+1)','cb');
end
