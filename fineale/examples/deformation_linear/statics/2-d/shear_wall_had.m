% H-adaptive analysis of stress distribution in a shear wall
function [utip]=shear_wall_had
    % Parameters:
    E=1e6;
    magn=1;
    graphics=true;
    Plot_errors=0;
    scale=10000;
    nu=1/3;
    convutip=23.97;
    a=1; b=5; c=1.5; d=4;
    bkgh=a;
    Region_definition={...
        ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(a) ' ' num2str(0)],...
        ['curve 2 line ' num2str(a) ' ' num2str(0) ' ' num2str(a) ' ' num2str(c)],...
        ['curve 3 line ' num2str(a) ' ' num2str(c) ' ' num2str(a+b) ' ' num2str(c)],...
        ['curve 4 line ' num2str(a+b) ' ' num2str(c) ' ' num2str(a+b) ' ' num2str(c+d)],...
        ['curve 5 line ' num2str(a+b) ' ' num2str(c+d) ' ' num2str(0) ' ' num2str(c+d)],...
        ['curve 6 line ' num2str(0) ' ' num2str(c+d) ' ' num2str(0) ' ' num2str(0)],...
        ['subregion 1  property 1 boundary  1  2 3 4 5 6']
        };
    Mesh_options =struct('quadratic',~true);
    Targetnel= 5000;
    convergence_rate=1.5;
    sigj=1;
                
    [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher(cat(2,Region_definition,{['m-ctl-point constant ' num2str(bkgh)]}), 1.0);
    
    for Adapt=1:3
        
        %     mesh{1}=fens;
        %     mesh{2}=fes;
        %     drawmesh(mesh); view(2); pause (1)
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
        ebc_fenids=fenode_select (fens,struct('box',[0,a,0,0],'inflate',a/10000));
        ebc_fixed=ones(1,length (ebc_fenids));
        ebc_comp=[];
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[a+b,a+b,0,c+d],'inflate',44/10000));
        ebc_fixed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fixed*0+1;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        % Number equations
        u   = numberdofs (u);
        % Assemble the system matrix
        K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
        % Load
        efemm = femm_deformation_linear (struct ('material',mater, ...
            'fes',subset(edge_fes,edge_groups{5}),...
            'integration_rule',gauss_rule (struct('dim',1,'order', 2))));
        fi=force_intensity(struct('magn',[0;-magn]));
        F =  distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        % Solve
        u = scatter_sysvec(u, K\F);
        
        
        nodal_stress = field_from_integration_points_spr(femm, geom, u, [], [], [], 'Cauchy', 1:3);
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
                %                             Plots stresses
                fld = field_from_integration_points_spr(femm, geom, u, [], [], [], 'Cauchy', sigj);
                nvals=fld.values;
                zeta=linspace(0,1,100)';
                p1=[0,1.001*c]; p2=p1+ [a+b,0];
                xi=ones(size(zeta,1),1)*p1+zeta*(p2-p1);
                fxi = simplex_grid_interpolation(geom.values,nvals,fes.conn,xi);
                plot(zeta,fxi,'-x'); hold on
                xlabel(' Normalized distance');
                ylabel([' Stress component \sigma_' num2str(sigj)]);
                drawnow
            end
            
        [fens,fes,groups,edge_fes,edge_groups] ...
            = targe2_mesher_adapt(Region_definition,fes.conn,geom.values,hests,1.0,Mesh_options);
        end
    end
