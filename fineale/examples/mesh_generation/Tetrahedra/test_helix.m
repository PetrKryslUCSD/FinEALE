% Mesh of a helix
[fens,fes] = T4_blocka(2.5*pi,0.95, 0.425, 50, 4, 2);
fens = transform_apply(fens,@(x, data) (x+ [0, 2.5, 0]), []);
climbPerRevolution= 2.3;
fens = transform_2_helix(fens,climbPerRevolution);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
