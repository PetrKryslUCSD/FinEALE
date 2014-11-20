% Mesh of two spherical layers.
function mksphelayrs
    R= 5;
    Ro=3*R;
    thickness =Ro-R;
    nlayers =13;
    [fens,fes]=H8_sphere(R,1);
    fes.label=1;
    V(fens,fes)
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    V(fens,fes)
    
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    
    function xyz= radially(xyz, layer)
        xyz= (R+layer/nlayers*thickness)*xyz/norm(xyz);
    end
    bg=mesh_boundary(fes);
    [fens1,fes1] = H8_extrude_Q4(fens,bg,nlayers,@radially);
    fes1.label=2;
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    
    bg=mesh_boundary(fes);
    %     l=gcell_select(fens,bg,struct ('facing', true, 'direction', [1,1,1]));
    function xyz= radially0(xyz, layer)
        xyz= (Ro+layer/nlayers*thickness)*xyz/norm(xyz);
    end
    [fens1,fes1] = H8_extrude_Q4(fens,bg,1,@radially0);
    fes1.label=0;
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, R/1000);
    fes=cat(fes1,fes2);
    
    l=fe_select(fens,fes,struct ('label', 1));
    gv=drawmesh({fens,mesh_boundary(subset(fes,l))},'fes','facecolor','red');
    l=fe_select(fens,fes,struct ('label', 2));
    gv=drawmesh({fens,mesh_boundary(subset(fes,l))},'gv',gv,'fes','facecolor','blue','facealpha',0.2);
    l=fe_select(fens,fes,struct ('label', 0));
    gv=drawmesh({fens,mesh_boundary(subset(fes,l))},'gv',gv,'fes','facecolor','y','facealpha',0.2);
    
    function V(fens,fes)
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        femm = femm_base (struct ('mater',[], 'fes',fes,...
            'integration_rule', gauss_rule (struct('dim',3,'order',4))));
        disp([' The volume is = '...
            num2str(integrate_function(femm,geom,inline('1'))) ])
    end
end