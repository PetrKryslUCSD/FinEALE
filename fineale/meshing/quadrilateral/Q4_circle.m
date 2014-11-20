function [fens,fes]=Q4_circle(radius,nrefine,thickness)
% Mesh of a quarter circle.
%
% function [fens,fes]=Q4_circle(radius,nrefine,thickness)
%
% Create a mesh of a quarter circle of "radius". The  mesh will consist of
% three quadrilateral elements if "nrefine==0", or more if "nrefine>0".
% "nrefine" is the number of bisections applied  to refine the mesh.
% The quadrilaterals have thickness "thickness".
%
% Examples: 
%     [fens,fes]=Q4_circle(1.35,1,1.0);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_sphere

    [fens,fes]=Q4_sphere(radius,nrefine,thickness);
    % apply transformation to project the locations of the nodes into the
    % plane x-y
    xyz=fens.xyz;
    for j=1:count(fens)
        xyz(j,1:2) = xyz(j,1:2)*((radius-norm(xyz(j,3)))+norm(xyz(j,3))/2)/radius;
    end
    fens.xyz=xyz(:,1:2);
end
