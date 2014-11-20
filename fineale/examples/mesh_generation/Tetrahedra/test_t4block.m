% Review of T4 block types
[fens,fes] = T4_blocka(3.0,1.2, 2, 4, 3, 2);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
title(strrep(('T4_blocka'),'_', '\_'))
pause(2); close all


[fens,fes] = T4_blockb(3.0,1.2, 2, 4, 3, 2);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
title(strrep(('T4_blockb'),'_', '\_'))
pause(2); close all

[fens,fes] = T4_blockca(3.0,1.2, 2, 4, 3, 2);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
title(strrep(('T4_blockca'),'_', '\_'))
pause(2); close all

[fens,fes] = T4_blockcb(3.0,1.2, 2, 4, 3, 2);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
title(strrep((' '),'_', '\_'))
pause(2); close all

