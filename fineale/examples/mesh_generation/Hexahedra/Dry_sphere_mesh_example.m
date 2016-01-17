function [fens,fes] = Dry_sphere_mesh
% This function creates the mesh of a sphere.
% The computational mesh is saved  in a Matlab .mat file, whose name is
% returned as output.
%
% The input parameters can be controlled  by selecting  values below.

    %     Now select the parameters:
    % Radius of the sphere
    R= 0.5;
    % Number of element edges through the radius of the sphere.
    %nradelems=13;
    nradelems=30;
    %     nradelems=10;
    
    [fens,fes]=H8_sphere_n(R,nradelems);
    fes.label=1;
    
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/nradelems/100);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/nradelems/100);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/nradelems/100);
    fes=cat(fes1,fes2);
    
    bfes= mesh_boundary(fes, []);
    disp(['Number of boundary cells = ' num2str(count(bfes))])
    %     drawmesh({fens,bfes},'fes','facecolor','yellow');
    %     drawmesh({fens,bfes},'fes','facecolor','none');
end

