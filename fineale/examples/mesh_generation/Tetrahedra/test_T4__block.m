% Conversion of T4 To T10 elements
X= [0,0,0;1,0,0;0,1,0;0,0,1];
fens=fenode_set(struct('xyz',X));
fes =fe_set_T4(struct ('conn',1:4));
[fens,fes] = T4_to_T10(fens,fes);

X=fens.xyz;
conn= fes.conn;
fes =fe_set_T4(struct ('conn',[7,6,3,10; 1,5,7,8; 5,2,6,9; 8,9,10,4; 10,7,6,9; 10,9,8,7; 5,6,7,9; 5,8,9,7]));
fens=fenode_set(struct('xyz',X(conn,:)));
% drawmesh({fens,fes},'fes','nodes','facecolor','red')
drawmesh({fens,fes},'nodes','fes','facecolor','red', 'shrink', 0.6)
% ix =gcell_select(fens,bg,...
%     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
% drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
axis off