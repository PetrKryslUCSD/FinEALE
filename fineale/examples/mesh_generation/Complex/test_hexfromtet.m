% Hexahedra from tetrahedra conversion.
[fens,fes] = T4_blocka(3.0,1.2, 2, 4, 6, 5);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')

tic;
[fens,fes] = T4_to_H8(fens,fes);
toc
% drawmesh({fens,fes},'fes','nodes','facecolor','red')
gv=drawmesh({fens,fes},'fes','facecolor','red');
interact(gv);
