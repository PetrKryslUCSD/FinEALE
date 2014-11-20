function [fens,fes]=Q4_circle_n(radius,nperradius,thickness)
% Mesh of a quarter circle with a given number of elements per radius.
%
% function [fens,fes]=Q4_circle_n(radius,nrefine,thickness)
%
% Create a mesh of 1/4 of the circle of "radius". The  mesh will consist of
% 3*(nperradius/2)^2 quadrilateral elements, which corresponds to 4*nperradius
% element edges per circumference. The parameter nperradius should be an even 
% number; if that isn't so is adjusted to by adding one. 
%
% Examples: 
%     [fens,fes]=Q4_circle_n(1.35,3,1.0);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_sphere_n

    [fens,fes]=Q4_sphere_n(radius,nperradius,thickness);
    % apply transformation to project the locations of the nodes into the
    % plane x-y
    xyz=fens.xyz;
    for j=1:count(fens)
        xyz(j,1:2) = xyz(j,1:2)*((radius-norm(xyz(j,3)))+norm(xyz(j,3))/2)/radius;
    end
    fens.xyz=xyz(:,1:2);
end
