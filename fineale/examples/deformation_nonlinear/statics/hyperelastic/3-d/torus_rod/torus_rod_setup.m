gv=graphic_viewer;
gv=reset(clear(gv,[]),[]);
set_graphics_defaults
set_label_defaults
cicl = gcell_select(fens, bgcells, struct('box', [0,0,-Inf,Inf,-Inf,Inf],'inflate',R/1000));
cbgcells =mesh_bdry(bgcells(cicl), struct('other_dimension',0.01));
cicl = gcell_select(fens, bgcells, struct('box', [-Inf,Inf,0,0,-Inf,Inf],'inflate',R/1000));
cbgcells =cat(2,cbgcells,mesh_bdry(bgcells(cicl), struct('other_dimension',0.01)));
cicl = gcell_select(fens, bgcells, struct('box', [-Inf,Inf,-Inf,Inf,0,0],'inflate',R/1000));
cbgcells =cat(2,cbgcells,mesh_bdry(bgcells(cicl), struct('other_dimension',0.01)));

options=struct ('x', geom, 'u', +0*u,'edgecolor','black','facecolor','black', 'shrink',1.0);
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end

xyz =gather (geom,(1:length(fens)),'values','noreshape');
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),xyz(:,2),xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),xyz(:,2),-xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),-xyz(:,2),-xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),-xyz(:,2),xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),xyz(:,2),xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),xyz(:,2),-xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),-xyz(:,2),-xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),-xyz(:,2),xyz(:,3)]));
for j=1:length(cbgcells),     draw(cbgcells(j),gv, options); end

options=struct ('x', geom, 'u', +0*u,'edgecolor','black','facecolor','white', 'shrink',1.0);
draw(sfeb,gv, options);
 gv=drawmesh( {fens,bgcells(icl)},'gv',gv,'facecolor','red');
camset(gv,[-122.4365 -155.3596  114.2994   -2.0019    1.5937    0.0791         0         0    1.0000 3.8487]);