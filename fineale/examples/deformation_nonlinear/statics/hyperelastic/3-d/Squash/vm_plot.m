set_graphics_defaults

gv=reset(clear(graphic_viewer,[]),[]);
% camset (gv,Cam);
cmap =cadcolors;
% draw(sfeb,gv, struct ('x', geom,'u',u, 'facecolor','red'));
% draw(sfeb,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
dT =u1;
Cauchy = [];
for j=1:6
    fld = field_from_integration_points(femm, geom, u1, [], 'Cauchy',j);
    nvals=fld.values;%min(nvals),max(nvals)
    Cauchy = [Cauchy,nvals];
end
nvals =sqrt(((Cauchy(:,1)-Cauchy(:,2)).^2+...
    (Cauchy(:,3)-Cauchy(:,2)).^2 +...
    (Cauchy(:,1)-Cauchy(:,3)).^2 +...
    6*(Cauchy(:,4).^2+Cauchy(:,5).^2 +Cauchy(:,6).^2))/2);
nvalsrange=[min(nvals),max(nvals)];
nvalsrange=[0,200];
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
draw(sfemm,gv, struct ('x', geom, 'u', +scale*u1,'colorfield',colorfield, 'edgecolor','black'));
% draw(feb,gv, struct ('x', geom, 'u', u,'facecolor','none',
% 'shrink',1.0));
% draw(efeb,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
colormap(cmap);
cbh=colorbar;
set(cbh,...
    'Position',[0.815 0.15 0.05 0.7],...
    'YLim',[0,1],...
    'YTick',[0,1],...
    'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
set(get(cbh,'XLabel'),'String','vm');
pause(0.5); camset(gv, [-3.7017   -7.6675    5.2828    0.5658    0.5304    0.3686         0         0    1.0000    7.9515]);
