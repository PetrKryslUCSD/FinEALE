% Test refinement of T4 tetrahedra
[fens,fes] = T4_cylinderdel(3.0,1.2, 2, 1);
[fens,fes] = T4_refine(fens,fes);
Volumes =t4util_volume (fes.conn,fens.xyz);
gv=drawmesh({fens,subset(fes,find(Volumes<=0))},'fes','facecolor','red','nodes')
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','gv',gv,'facecolor','none')
% drawmesh({fens,bg},'fes','facecolor','red')
% % ix =fe_select(fens,bg,...
% %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
% % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
