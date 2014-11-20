% Test generation of a hollow cylinder.
[fens,fes] = H8_hollow_cylinder(0.5, 0.2, 0.7, 14, 2, 5)
bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','facecolor','y');

