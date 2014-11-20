% Test conversion of T4 tetrahedra into hexahedra and quadratic tetrahedra.

[fens,fes] = T4_blocka(3.0,1.2, 2, 4, 3, 2);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
title(strrep(('T4_blocka'),'_', '\_'))


[fens,fes] = T4_to_H8(fens,fes);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','y')
title('Converted from T4 to H8')

figure
[fens,fes] = T4_blockb(3.0,1.2, 2, 4, 3, 2);
[fens,fes] = T4_to_T10(fens,fes);
% [fens,fes] = T10_to_H8(fens,fes);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','b')
title('Converted from T4 to T10')



