% Make quadrilateral surface  mesh with given subdivision.
function [xyz,conn]=make_spherical_surface_n(Radius,nperradius)
    if (~ exist ('Radius') )
        Radius=300;
    end
    if (~ exist ('nperradius') )
        nperradius=4;
    end
    tol=Radius/nperradius/2/1000;
    [fens,fes]=Q4_sphere_n(Radius,nperradius,1);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    drawmesh({fens,fes},'fes','facecolor','red')
    xyz=fens.xyz;
    conn =fes.conn;
end