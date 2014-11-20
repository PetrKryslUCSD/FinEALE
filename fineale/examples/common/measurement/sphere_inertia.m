% Compute the  geometry characteristics of a sphere.
function sphere_inertia(Radius,nsubdivisions)
    if (~ exist ('Radius') )
        Radius=3;
    end
    if (~ exist ('nsubdivisions') )
        nsubdivisions=2;
    end
    R=Radius;
    tol=R/nsubdivisions/2/1000;
    [fens,fes]=H8_sphere(R,nsubdivisions);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
    fes=cat(fes1,fes2);
    CG1= 10*[1,-2,3];
    fens = transform_apply(fens,@(x,d)x+CG1, []);
    drawmesh({fens,fes},'fes','facecolor','red');
     
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    femm = femm_base (struct ('mater',[], 'fes',fes,...
        'integration_rule', gauss_rule (struct('dim',3,'order',2))));
    V=integrate_function(femm,geom,@(x)(1));
    disp( [' Volume = '  num2str(V) ' (' num2str(4/3*pi*R^3) ')'])
    S=integrate_function(femm,geom,@(x)(x));
    CG=S/V;
    disp( [' Center of gravity = '  num2str(CG) ' (' num2str(CG1) ')'])
    I=integrate_function(femm,geom,@(x)(norm(x-CG)^2*eye(3)-(x-CG)'*(x-CG)));
    disp( [' Tensor of inertia = ' ])
    I
    2/5*V*R^2*eye(3)
end