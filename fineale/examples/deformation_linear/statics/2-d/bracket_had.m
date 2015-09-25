% H-adaptive analysis of stress distribution in a shear wall
function [utip]=shear_wall_had
% Parameters:
E=70e9;
magn=1e7;
graphics=true;
Plot_errors=0;
scale=100;
nu=0.33;
a=0.08; b=0.12; c=0.04; d=0.06; e= 0.05; R1= 0.01; R2= 0.005; f= 0.02; t=  0.005;
bkgh=R2;
Region_definition={...
    ['curve 1 line ' num2str(c) ' ' num2str(0) ' ' num2str(a) ' ' num2str(0)],...
    ['curve 2 line ' num2str(a) ' ' num2str(0) ' ' num2str(a) ' ' num2str(b)],...
    ['curve 3 line ' num2str(a) ' ' num2str(b) ' ' num2str(0) ' ' num2str(b)],...
    ['curve 4 line ' num2str(0) ' ' num2str(b) ' ' num2str(0) ' ' num2str(d)],...
    ['curve 5 line ' num2str(0) ' ' num2str(d) ' ' num2str(c) ' ' num2str(d)],...
    ['curve 6 arc ' num2str(c) ' ' num2str(d) ' ' num2str(c) ' ' num2str(d-2*R1) ' center ' num2str(0.99*c) ' ' num2str(d-R1)],...
    ['curve 7 line ' num2str(c) ' ' num2str(d-2*R1) ' ' num2str(c) ' ' num2str(0)],...
    ['curve 8 Circle Center ' num2str(f) ' ' num2str(b-f) ' Radius ' num2str(R2)],...
    ['curve 9 Circle Center ' num2str(a-f) ' ' num2str(b-f) ' Radius ' num2str(R2)],...
    ['subregion 1  property 1 boundary  1  2 3 4 5 6 7 hole -8 -9 ']
    };
Mesh_options =struct('quadratic',~true);
Targetnel= 7000;
convergence_rate=1.5;
sigj=3;

[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher(cat(2,Region_definition,{['m-ctl-point constant ' num2str(bkgh)]}), t);

for Adapt=1:3
    
    %     mesh{1}=fens;
    %     mesh{2}=fes;
    %     drawmesh(mesh); view(2); pause (1)
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
    ebc_fenids=connected_nodes (subset(edge_fes,cat(2,edge_groups{8:9})));
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
    efemm = femm_deformation_linear (struct ('material',mater, ...
        'fes',subset(edge_fes,edge_groups{7}),...
        'integration_rule',gauss_rule (struct('dim',1,'order', 2))));
    fi=force_intensity(struct('magn',[magn;0]));
    F =  distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
    % Solve
    u = scatter_sysvec(u, K\F);
    
    
    nodal_stress = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', 1:3);
    elerrs = flux_L2_error (femm, geom, u, [], nodal_stress);
    total_err=sqrt(sum(elerrs.^2));
    targeterr=sqrt(total_err^2/Targetnel);
    
    [hcurs, hests] =T3_mesh_sizes(fes.conn,geom.values,targeterr,elerrs,convergence_rate);
    
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
            gv=reset (gv,struct('limits', [-0.1*a,1.1*a,-0.1*b,1.1*b]));
            cmap=jet;
            fld = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', sigj);
            nvals=fld.values;
            dcm=data_colormap(struct ('range',[min(nvals),max(nvals)], 'colormap',cmap));%[min(nvals),max(nvals)]
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
            draw(femm, gv, struct ('x',geom,'u', scale*u, 'colorfield',colorfield));
            view (2)
            lighting  none;
            nvalsrange=[min(nvals),max(nvals)];
            cbh=colorbar;
            set(cbh,...
                'Position',[0.8 0.15 0.05 0.7],...
                'YLim',[0,1],...
                'YTick',[0,1],...
                'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
            set(get(cbh,'XLabel'),'String',['\sigma_' num2str(sigj)]);
            pause(0.1);
        end
        
        [fens,fes,groups,edge_fes,edge_groups] ...
            = targe2_mesher_adapt(Region_definition,fes.conn,geom.values,hests,1.0,Mesh_options);
    end
end
