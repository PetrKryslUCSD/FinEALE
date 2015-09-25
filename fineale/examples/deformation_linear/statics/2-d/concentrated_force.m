% H-adaptive analysis of stress distribution in a shear wall
function [utip]=shear_wall_had
% Parameters:
E=70e9;
magn=1e3;
graphics=true;
Plot_errors=0;
scale=200;
nu=0.33;
a=0.08; b=0.12;  t=  0.005;
bkgh=a/15;
Region_definition={...
    ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(a) ' ' num2str(0)],...
    ['curve 2 line ' num2str(a) ' ' num2str(0) ' ' num2str(a) ' ' num2str(b)],...
    ['curve 3 line ' num2str(a) ' ' num2str(b) ' ' num2str(0) ' ' num2str(b)],...
    ['curve 4 line ' num2str(0) ' ' num2str(b) ' ' num2str(0) ' ' num2str(0)],...
    ['subregion 1  property 1 boundary  1  2 3 4 '],...
  ['m-ctl-point 1 xy ' num2str(a) ' ' num2str(b) ' near ' num2str(bkgh/2000) ' influence ' num2str(bkgh/100)],...
    ['m-ctl-point constant ' num2str(bkgh)],...
    };
% ['m-ctl-point 1 xy ' num2str(a) ' ' num2str(b) ' near ' num2str(bkgh/2000) ' influence ' num2str(bkgh/100)],...
    

[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher(Region_definition, t);


%     mesh{1}=fens;
%     mesh{2}=fes;
%     drawmesh(mesh,'fes'); view(2); pause (1)
disp(['Number of elements: ',num2str(count(fes))])

%     Material
% Material
prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_linear_biax (struct('property',prop, ...
    'reduction','stress'));
% Finite element block
femm = femm_deformation_linear(struct ('material',mater, 'fes',fes,...
    'integration_rule',tri_rule (struct('npts',1))));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero outfunction five
% Apply EBC's
ebc_fenids=connected_nodes (subset(edge_fes,cat(2,edge_groups{1})));
ebc_fixed=ones(1,length (ebc_fenids));
ebc_comp=[];
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);
% Assemble the system matrix
K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
% Load
cn=fenode_select(fens,struct('box', bounding_box ([a,b]), 'inflate',a/1e5));
efemm = femm_deformation_linear (struct ('material',mater, ...
    'fes',fe_set_P1(struct('conn',cn)),...
    'integration_rule',point_rule ));
fi=force_intensity(struct('magn',[magn;0]));
F =  distrib_loads(efemm, sysvec_assembler, geom, u, fi, 3);
norm(F)
% Solve
u = scatter_sysvec(u, K\F);
cu=gather_values(u,cn);
disp(['Displacement under the force =', num2str(cu(1))])

% Plot
if graphics
    %             Plot errors
    if (Plot_errors)
        gv=graphic_viewer;
        gv=reset (gv,struct('limits', [-a,a+b,0,c+d+a]));
        cmap=cadcolors;
        nvals=sqrt(elerrs);
        dcm=data_colormap(struct ('range',[min(nvals),max(nvals)], 'colormap',cmap));%[min(nvals),max(nvals)]
        for i=1:count(fes)
            color = map_data (dcm, nvals(i));
            draw(subset(fes,i), gv, struct ('x',geom,'u', u, 'facecolor',color));
        end
        view (2)
        lighting  none;
        nvalsrange=[min(nvals),max(nvals)];
        cbh=colorbar;
        set(cbh,...
            'Position',[0.8 0.15 0.05 0.7],...
            'YLim',[0,1],...
            'YTick',[0,1],...
            'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
        set(get(cbh,'XLabel'),'String',['Error']);
    else
        %                 %                             Plots stresses
        %                 fld = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', sigj);
        %                 nvals=fld.values;
        %                 zeta=linspace(0,1,100)';
        %                 p1=[0,1.001*c]; p2=p1+ [a+b,0];
        %                 xi=ones(size(zeta,1),1)*p1+zeta*(p2-p1);
        %                 fxi = simplex_grid_interpolation(geom.values,nvals,fes.conn,xi);
        %                 plot(zeta,fxi,'-x'); hold on
        %                 xlabel(' Normalized distance');
        %                 ylabel([' Stress component \sigma_' num2str(sigj)]);
        %                 drawnow
        gv=graphic_viewer;
        gv=reset (gv,struct('limits', [-0.1*a,1.63*a,-0.1*b,1.1*b]));
        cmap=jet;
        %             fld = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', sigj);
        %             nvals=fld.values;
        nvals=sqrt(u.values(:,1).^2+u.values(:,2).^2);
        dcm=data_colormap(struct ('range',[min(nvals),max(nvals)], 'colormap',cmap));%[min(nvals),max(nvals)]
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
        draw(femm, gv, struct ('x',geom,'u', scale*u, 'colorfield',colorfield));
        view (2)
        lighting  none;
        nvalsrange=[min(nvals),max(nvals)];
        cbh=colorbar;
        set(cbh,...
            'Colormap',cmap,...
            'Position',[0.8 0.15 0.05 0.7],...
            'YLim',[0,1],...
            'YTick',[0,1],...
            'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
        set(get(cbh,'XLabel'),'String',['|u|']);
        pause(0.1);
    end
    
    
end
