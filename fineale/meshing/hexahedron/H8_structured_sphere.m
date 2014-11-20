function [fens,fes] = H8_structured_sphere(R, Mesh_size)
% Mesh of an entire sphere.
%
% function [fens,fes] = H8_structured_sphere(R, Mesh_size)
%
% R= radius of the sphere, 
% Mesh_size= mesh size, 
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%    R= 0.5; Mesh_size=R/4;
%    [fens,fes] = H8_structured_sphere(R, Mesh_size)
%    drawmesh ({fens,fes}, 'fes', 'facecolor','y');
%
% See also: H8_sphere_n


   if (~exist('R','var'))
        R= 0.5;
    end
    if (~exist('Mesh_size','var'))
        Mesh_size=R/2;
    end
    
    nperradius=round(R/ Mesh_size);
    if (nperradius<1)
    	nperradius=1;
    end
    tol=R/nperradius/1000;% geometrical tolerance
    
    % Mesh
    [fens,fes]=H8_sphere_n(R,nperradius);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
end