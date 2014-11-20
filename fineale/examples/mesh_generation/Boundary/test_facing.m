% Test selection of elements based on which way they are facing.
[fens,fes] = T4_blocka(0.5*pi,0.95, 0.425, 5, 4, 2);

fens = transform_apply(fens,@(x, data) (x+ [0, 2.5, 0]), []);
climbPerRevolution= 0;
fens = transform_2_helix(fens,climbPerRevolution);
bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','facecolor','none')


l=fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,1]));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y')
l=fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,-1]));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','r')
l=fe_select(fens,bg,struct ('facing', true, 'direction', @(x)(-x),'tolerance',eps));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','b')
l=fe_select(fens,bg,struct ('facing', true, 'direction', @(x)(+[x(1:2),0]),'tolerance',0.01));
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','c')
