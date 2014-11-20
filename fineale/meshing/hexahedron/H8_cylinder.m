function [fens,fes] = H8_cylinder(Radius, Length, nref, nL)
% H8 mesh of a solid  cylinder with a given edges per radius.
% 
% function  [fens,fes] = H8_cylinder(Radius, Length, nref, nL)
%   
% Make H8 mesh of a solid  cylinder with a given number of refinements 
% by bisection of the elements within the cross-section.
%
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% Examples: 
%     [fens,fes] = H8_cylinder(1.0, 3.0, 1, 3);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_cylinder_n
%
     tol=min( [Length/nL,Radius/2/nref] )/100;
    [fens,fes] = Q4_circle(Radius, nref, 1.0);
    fens2 = transform_apply(fens,@(x,d)([-x(2),x(1)]),  0);
    [fens,fes,fes2] = merge_meshes(fens, fes, fens2, fes, tol);
    fes = cat(fes,fes2);
    fens2 = transform_apply(fens,@(x,d)([-x(1),-x(2)]),  0);
    [fens,fes,fes2] = merge_meshes(fens, fes, fens2, fes, tol);
    fes = cat(fes,fes2);
    [fens,fes] = H8_extrude_Q4(fens,fes,nL,@(x,k)([x(1),x(2),k*Length/nL]));
end