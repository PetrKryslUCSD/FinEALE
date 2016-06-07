% Test generation of a cylinder mesh.
 [fens,fes] = H8_quarter_cylinder_n(30, 80, 4, 7);
count(fens)
bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','facecolor','none');
