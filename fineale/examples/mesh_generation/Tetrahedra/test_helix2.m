% Mesh of circular cross-section helix.
[fens,fes] = T4_cylinderdel(3.0,1.2, 8, 5);
fens = transform_apply(fens,@(x, data) (x+ [0, -2.5, 0]), []);
climbPerRevolution= 2.3;
fens = transform_2_helix(fens,climbPerRevolution);
bg=mesh_boundary(fes);
drawmesh({fens,bg},'fes','facecolor','red')
% ix =gcell_select(fens,bg,...
%     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
% drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
