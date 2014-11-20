% Test selection of elements based on which way they are facing.
[fens,fes] = T4_blockca(0.5*pi,0.95, 0.425, 4, 2, 3);
fens = transform_apply(fens,@(x, data) (x+ [0, 2.5, 0]), []);
climbPerRevolution= 4.3;
fens = transform_2_helix(fens,climbPerRevolution);

bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','nodes','facecolor','none')


l=fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,1]));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y')
l=fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,-1]));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','r')
