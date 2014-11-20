% Make one octant of a sphere with core and mantle.
function mkoctsphelayrs
    R= 5;
    Ro=3*R;
    thickness =Ro-R;
    nlayers =13;
    [fens,fes]=H8_sphere(R,1);
    fes.label=1;
    function xyz= radially(xyz, layer)
        xyz= (R+layer/nlayers*thickness)*xyz/norm(xyz);
    end
    bg=mesh_boundary(fes);
    l=fe_select(fens,bg,struct ('facing', true, 'direction', [1,1,1]));
    [fens1,fes1] = H8_extrude_Q4(fens,subset(bg,l),nlayers,@radially);
    fes1.label=2;
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    l=fe_select(fens,fes,struct ('label', 1));
    gv=drawmesh({fens,subset(fes,l)},'fes','facecolor','red')
    l=fe_select(fens,fes,struct ('label', 2));
    gv=drawmesh({fens,subset(fes,l)},'gv',gv,'fes','facecolor','blue')
end