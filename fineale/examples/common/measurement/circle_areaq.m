% Example of measurement of geometry characteristics of planar figures, quadratic triangles.
% In this example we are measuring the area of a circle using its triangulation.
% The triangles are quadratic, curved so that the circular boundary is 
% approximated by a parabolic arc.

R=25;
center =[13 4.2];
exact_area =pi*R^2;
errors = []
mesh_sizes = [32, 16, 8, 4];
for mesh_size = mesh_sizes
    [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
        ['curve 1 circle  Center ' num2str(center) '  radius ' num2str(R)],...
        'subregion 1  property 1 boundary 1',...
        ['m-ctl-point constant ' num2str(mesh_size)]
        }, 1.0, struct('quadratic',true));
    el =edge_groups{1};
    s=subset(edge_fes,el);
    conn =s.conn;
    xyz=fens.xyz;
    for i=1:size(conn,1)
        r =xyz(conn(i,3),:)-center;
        xyz(conn(i,3),:) =R*r/norm(r)+center;
    end 
    fens.xyz=xyz;
    gv =drawmesh({fens,fes},'fes','shrink', 1);
    view(2)
    geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
    femm = femm_base (struct ('mater',[], 'fes',fes,...
        'integration_rule', tri_rule(struct('npts',6))));
    approximate_area =integrate_function(femm,geom,inline('1','x'));
    errors = [errors abs(exact_area-approximate_area)]
    disp([' The area is = ' num2str(approximate_area) ', to be compared with ' num2str(exact_area)])
    pause(1);   clear (gv); close(gcf);
end
loglog(mesh_sizes, errors, 'bo-'); grid on
assignin('caller','fineale_test_passed',((norm([1.963492595470922e+003]-approximate_area))<1e-9))