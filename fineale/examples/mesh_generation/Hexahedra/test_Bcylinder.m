% Test generation of a cylinder mesh.
 [fens,fes] = H8_cylinder(0.195, 0.5, 1, 3);
count(fens)
bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','facecolor','none');


% l=fe_select(fens,bg,struct ('facing', true, 'direction', @(x)(-[x(1:2),0]),'tolerance',eps));
% gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','b')


l=[fe_select(fens,bg,struct ('facing', true, 'direction', @(x)(-[x(1:2),0]),'tolerance',eps)),...
    fe_select(fens,bg,struct ('facing', true, 'direction', @(x)([x(1:2),0]),'tolerance',eps))];
gv =drawmesh({fens,subset(bg,setdiff(1:count(bg),l))},'gv',gv,'fes','facecolor','b');
l=[fe_select(fens,bg,struct ('facing', true, 'direction', @(x)(-[x(1:2),0]),'tolerance',eps))];
gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y');
