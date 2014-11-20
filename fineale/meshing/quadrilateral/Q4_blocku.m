function [fens,fes] = Q4_blocku(Length,Width,nL,nW,options)
% Unstructured mesh of a rectangle, Q4 elements.
%
% function [fens,fes] = Q4_blocku(Length,Width,nL,nW,options)
%
% Range =<0,Length> x <0,Width>
% Divided into Quadrilateral (Q4) elements: nL, nW in the first, second (x,y).
% This function  produces an unstructured mesh of Q4 quadrilaterals by
% generating a triangular mesh and then using trisection to split the
% triangles into quadrilaterals.
%
% Examples: 
%     [fens,fes] = Q4_blocku(3.5,1.75,2,3,struct('other_dimension',1))
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: T3_block2du, T3_to_Q4, Q4_refine


    if ~isstruct(options)
        other_dimension = options; clear options;
        options.other_dimension = other_dimension;
    end
    [fens,fes] = T3_blocku(Length,Width,1,1,options);
    [fens,fes] = T3_to_Q4(fens,fes,options);
    n=2;
    while (2*n < max([nL,nW]))
        [fens,fes] = Q4_refine(fens,fes);
        n =n*2;
    end
end
