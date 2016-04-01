set_graphics_defaults
Cam =[-0.8872   -1.3668    1.0361    0.2000    0.0500    0.0050         0         0    1.0000    7.0899 ];
gv=reset(clear(graphic_viewer,[]),[]);
camset (gv,Cam);
cmap =cadcolors;
% draw(sfeb,gv, struct ('x', geom,'u',u, 'facecolor','red'));
% draw(sfeb,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
dT =u1;
% Cauchy = [];
% for j=1:6
%     fld = field_from_integration_points(femm, geom, u1, [], 'Cauchy',j);
%     nvals=fld.values;%min(nvals),max(nvals)
%     Cauchy = [Cauchy,nvals];
% end
% nvals =sqrt(((Cauchy(:,1)-Cauchy(:,2)).^2+...
%     (Cauchy(:,3)-Cauchy(:,2)).^2 +...
%     (Cauchy(:,1)-Cauchy(:,3)).^2 +...
%     6*(Cauchy(:,4).^2+Cauchy(:,5).^2 +Cauchy(:,6).^2))/2);
% nvalsrange=[min(nvals),max(nvals)];
nvalsrange=[0,1e5];

dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
% colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
% draw(sfemm,gv, struct ('x', geom, 'u', +scale*u1,'colorfield',colorfield, 'edgecolor','black'));
boundaryfes = mesh_boundary (femm.fes,[]);
    oon_boundary =zeros(count(fens),1);
    oon_boundary(boundaryfes.conn(:)) =1;
    fakefemm =femm;
    for feix=1:count(femm.fes)
        fakefemm.fes =subset(femm.fes,feix);
        if (sum(oon_boundary(fakefemm.fes.conn(:)))>0)
            % Create the color field
            Cauchy = [];
            for j=1:6
                fld = field_from_integration_points(fakefemm, geom, u, [], 'Cauchy',j);
                nvals=fld.values;
                Cauchy = [Cauchy,nvals];
            end
            nvals =sqrt(((Cauchy(:,1)-Cauchy(:,2)).^2+...
                (Cauchy(:,3)-Cauchy(:,2)).^2 +...
                (Cauchy(:,1)-Cauchy(:,3)).^2 +...
                6*(Cauchy(:,4).^2+Cauchy(:,5).^2 +Cauchy(:,6).^2))/2);
            nvalsrange=[min([nvalsrange,min(nvals)]),max([nvalsrange,max(nvals)]),];
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
                map_data(dcm, nvals)));
            draw(fakefemm, gv, struct ('x',geom, 'u',u,...
                'colorfield',colorfield, 'shrink',1.0));
        end
    end
    nvalsrange
nvalsrange=[0,1e5];
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
pause(0.5); 
