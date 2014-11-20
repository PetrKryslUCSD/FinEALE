% Example of measurement of geometry characteristics of planar figures.
% In this example we are measuring the area of a circle using its triangulation.
exact_area =pi*25^2;
errors = []
mesh_sizes = [32, 16, 8, 4, 2];
for mesh_size = mesh_sizes
    [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
        'curve 1 circle  Center 13 4.2 radius 25',...
        'subregion 1  property 1 boundary 1',...
        ['m-ctl-point constant ' num2str(mesh_size)]
        }, 1.0);
    gv=drawmesh({fens,fes},'fes','shrink', 1.0);
    view(2)
    geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
    femm = femm_base (struct ('mater',[], 'fes',fes,...
        'integration_rule', tri_rule(struct('npts',1))));
    approximate_area =integrate_function(femm,geom,@(x)1);% if anonymous functions are not available, use inline('1','x')
    errors = [errors abs(exact_area-approximate_area)]
    disp([' The area is = ' num2str(approximate_area) ', to be compared with ' num2str(exact_area)])
    pause(1);   clear (gv); close(gcf);
end
loglog(mesh_sizes, errors, 'bo-'); grid on