set_graphics_defaults
set_label_defaults
gv=graphic_viewer;
gv=reset(clear(gv,[]),[]);
cmap =cadcolors;
%             draw(sfeb,gv, struct ('x', geom,'u',u, 'facecolor','red'));
%             draw(sfeb,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
fld = field_from_integration_points(feb, geom, u1, [], 'Cauchy',3);
nvals=get(fld,'values');%min(nvals),max(nvals)
nvalsrange=0.5*[min(nvals),max(nvals)];
nvalsrange=[-0.5, 0.5]*1e6;
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
xyz =gather (geom,(1:length(fens)),'values','noreshape')+...
    gather (u,(1:length(fens)),'values','noreshape');
options=struct ('x', geom, 'u', +0*u,'colorfield',colorfield,'edgecolor','black', 'shrink',1.0);
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),xyz(:,2),xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),xyz(:,2),-xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),-xyz(:,2),-xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [xyz(:,1),-xyz(:,2),xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),xyz(:,2),xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),xyz(:,2),-xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),-xyz(:,2),-xyz(:,3)]));
draw(sfeb,gv, options);
options.x=field(struct ('name', ['mgeom'], 'data', [-xyz(:,1),-xyz(:,2),xyz(:,3)]));
draw(sfeb,gv, options);
% draw(sfeb,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
% draw(efeb,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
camset(gv,[-122.4365 -155.3596  114.2994   -2.0019    1.5937    0.0791         0         0    1.0000 3.8487]);
colormap(cmap);
cbh=colorbar;
set(cbh,...
    'Position',[0.0382    0.1381    0.0500    0.7000],...
    'YLim',[0,1],...
    'YTick',[0,1],...
    'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
set(get(cbh,'XLabel'),'String','\sigma_z');