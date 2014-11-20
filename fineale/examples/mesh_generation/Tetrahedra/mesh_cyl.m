% Delaunay mesh of a cylinder with T4 and T10 tetrahedra
R= 1.2; L= 3.0;
[fens,fes] = T4_cylinderdel(L,R, 2, 1);
bg=mesh_boundary(fes);
gv=drawmesh({fens,bg},'fes','facecolor','y')
axis on
labels('X', 'Y', 'Z')
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
show_field_as_markers(gv, struct('x',geom,'u',0*geom,'marker','o','markersize',4 ))
[fens,fes] = T4_to_T10(fens,fes);
bg=mesh_boundary(fes);
l=[fe_select(fens,bg,struct ('facing', true, 'direction', @(x)([0,x(2:3)]),'tolerance',eps))];
% gv =drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y','facealpha',0.8);
axis on
labels('X', 'Y', 'Z')
pause (1)

xyz=fens.xyz;
sbg  =subset(bg,l);
conn=sbg.conn;

for je=1:size(conn,1 )
    for jk=1:size(conn,2)
        x1=xyz(conn(je,jk),:);
        xyz(conn(je,jk),:)= [x1(1),x1(2)/norm(x1([2,3]))*R,x1(3)/norm(x1([2,3]))*R];
    end
end
fens.xyz=xyz;
gv=drawmesh({fens,bg},'fes','facecolor','g')
axis on
labels('X', 'Y', 'Z')
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
show_field_as_markers(gv, struct('x',geom,'u',0*geom,'marker','o','markersize',4 ))