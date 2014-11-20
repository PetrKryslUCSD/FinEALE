function [fens,fes]=Q4_sphere(radius,nrefine,thickness)
% Surface mesh of 1/8 of a sphere. 
%
% function [fens,fes]=Q4_sphere(radius,nrefine,thickness)
%
% Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
% three quadrilateral elements if "nrefine==0", or more if "nrefine>0".
% "nrefine" is the number of bisections applied  to refine the mesh.
% The quadrilaterals have thickness "thickness".
%
% Examples: 
% [fens,fes]=Q4_sphere(77.1,2,1.0);
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_sphere_n

    fens =fenode_set(struct('xyz',[[1, 0, 0];[0, 1, 0];[0, 0, 1];[1, 1, 0];[0, 1, 1];[1, 0, 1];[1, 1, 1]]));
    fes=fe_set_Q4(struct ('conn',[[1, 4, 7, 6];[4, 2, 5, 7]; [3, 6, 7, 5]],'other_dimension', thickness));
    fens= onto_sphere(fens, radius);
    for i=1:nrefine
        [fens,fes]=Q4_refine(fens,fes);
        fens= onto_sphere(fens, radius);
    end
end
