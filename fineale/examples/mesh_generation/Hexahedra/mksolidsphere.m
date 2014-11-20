% Mesh of one octant of a solid sphere.
function mksolidsphere
    R= 300;
    Ri= 1;
    thickness =R-Ri;
    nlayers =13;
    [fens,fes]=H8_sphere(R,1);
    drawmesh({fens,fes},'fes','facecolor','red')
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    femm = femm_base (struct ('mater',[], 'fes',fes,...
        'integration_rule', gauss_rule (struct('dim',3,'order',4))));
    disp([' The volume is = '...
        num2str(integrate_function(femm,geom,inline('1'))) ...
        ', to be compared with ' num2str(4/3*pi*R^3/8)])
end