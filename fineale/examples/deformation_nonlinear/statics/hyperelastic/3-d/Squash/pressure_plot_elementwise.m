gv=reset(clear(graphic_viewer,[]),[]);
camset (gv,Cam);
% draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
% draw(sfemm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
% fld = field_from_integration_points(femm, geom, u1, [], 'pressure',1);
% nvals=fld.values;%min(nvals),max(nvals)
nvalsrange=[0,320];
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
% colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
% draw(sfemm,gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink',1.0));
% draw(femm,gv, struct ('x', geom, 'u', u,'facecolor','none', 'shrink',1.0));
% draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
boundaryfes = mesh_boundary (femm.fes,[]);
    oon_boundary =zeros(count(fens),1);
    oon_boundary(boundaryfes.conn(:)) =1;
    fakefemm =femm;
    for feix=1:count(femm.fes)
        fakefemm.fes =subset(femm.fes,feix);
        if (sum(oon_boundary(fakefemm.fes.conn(:)))>0)
            % Create the color field
            Cauchy = [];
            fld = field_from_integration_points(fakefemm, geom, u, [], 'pressure',1);
                nvals=fld.values;
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
                map_data(dcm, nvals)));
            draw(fakefemm, gv, struct ('x',geom, 'u',u,...
                'colorfield',colorfield, 'shrink',1.0));
        end
    end
colormap(cmap);
cbh=colorbar;
set(cbh,...
    'Position',[0.815 0.15 0.05 0.7],...
    'YLim',[0,1],...
    'YTick',[0,1],...
    'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
set(get(cbh,'XLabel'),'String','pressure');
pause(0.5); Cam =camget(gv);
axis off
set_pub_defaults