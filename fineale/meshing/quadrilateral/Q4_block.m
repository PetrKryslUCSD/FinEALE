function [fens,fes] = Q4_block(Length,Width,nL,nW,options)
% Mesh of a rectangle, Q4 elements
%
% function [fens,fes] = Q4_block(Length,Width,nL,nW,options)
%
% Rectangle <0,Length> x <0,Width>
% Divided into elements: nL, nW in the first, second (x,y).
% options = structure with fields recognized by the constructor of the
%       Q4 finite element
% 
% Examples: 
% [fens,fes] = Q4_block(3.5,1.75,2,3,struct('other_dimension',1))
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%  
% See also: Q4_blockx, fe_set_Q4

    [fens,fes] = Q4_blockx(linspace(0,Length,nL+1)',linspace(0,Width,nW+1)',options);
end
