function [fens,fes] = T4_sphere(R,nr)
% Tetrahedral (t4) Delaunay mesh of a sphere.
%
% function [fens,fes] =  T4_sphere(R,nr)
%
% R= radius, 
% nr= number of element edges radially, call for H8_sphere_n
%  
% Examples: 
% [fens,fes] = T4_sphere(3.1,4);
% figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
% See also: H8_sphere_n

    tol=R/nr/100;
    [fens,fes]=H8_sphere_n(R,nr);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    T = delaunayn(fens.xyz);
    fes=fe_set_T4(struct ('conn',T(:,:)));
end

%compute volume of a tetrahedron
% Given the 4x3 vertex coordinate matrix V of a tetrahedron, TETVOL(V)
% returns the volume of the tetrahedron.
function vol = tetvol(v)
    vol = det([v(2,:)-v(1,:);v(3,:)-v(1,:);v(4,:)-v(1,:)])/6;
    if abs (vol) < 0.1
        warning (' sliver?')
    end
    return;
end
