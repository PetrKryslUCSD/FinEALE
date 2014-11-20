function [fens,fes] = H8_block(Length,Width,Height,nL,nW,nH)
% Mesh of a 3-D block of H8 finite elements
%
% function [fens,fes] = H8_block(Length,Width,Height,nL,nW,nH)
%
% Arguments:
% Length,Width,Height= dimensions of the mesh in Cartesian coordinate axes,
% smallest coordinate in all three directions is  0 (origin)
% nL,nW,nH=number of elements in the three directions
%
% Range in xyz =<0,Length> x <0,Width> x <0,Height>
% Divided into elements: nL, nW, nH in the first, second, and
% third direction (x,y,z). Finite elements of type H8.
%
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% Examples: 
%     [fens,fes] = H8_block(2,3,4, 2,6,4);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_block_u, H8_blockx, H8_composite_plate

    [fens,fes] = H8_blockx(linspace(0,Length,nL+1)',linspace(0,Width,nW+1)',linspace(0,Height,nH+1)');
end
