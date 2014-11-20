function [fens,fes] = T4_blockb(Length,Width,Height,nL,nW,nH)
    % Tetrahedral mesh of a rectangular block; orientation 'b'.
    %
    % function [fens,fes] = T4_blocka(Length,Width,Height,nL,nW,nH)
    %
    % The mesh is produced by splitting each logical  rectangular cell into six
    % tetrahedra.
    % Range =<0,Length> x <0,Width> x <0,Height>
    % Divided into elements: nL, nW, nH in the first, second, and
    % third direction (x,y,z).
    %
    % Output:
    % fens= finite element node set
    % fes = finite element set
    %
    % Examples:
    %         [fens,fes] = T4_blockb(2,3,4, 2,6,4);
    %         figure; drawmesh({fens,fes},'fes','facecolor','g'); hold on
    %
    % See also: T4_blockx,  T4_blocka,  T4_blockca,  T4_blockcb
    [fens,fes] = T4_blockx(linspace(0,Length,nL+1)',linspace(0,Width,nW+1)',linspace(0,Height,nH+1)','b');
end
