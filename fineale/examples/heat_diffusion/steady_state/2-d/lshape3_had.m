% Adaptive solution on L-shaped region.
% L-shaped region (one quarter of a rectangular region with a rectangular
% hole).   Uniform heat source, temperature prescribed on the boundary.
%
%
function lshape3_had
    kappa=[0.2 0; 0 0.2]; % conductivity matrix
    Q=0.013; % uniform heat source
    
    region_definition={'REGEXT 0 0 48 48',...
        ['curve 1 line 20 0 48 0'],...
        ['curve 2 line 48 0 48 48'],...
        ['curve 3 line 48 48 0 48'],...
        ['curve 4 line 0 48 0 13'],...
        ['curve 5 line 0 13 20 13'],...
        ['curve 6 line 20 13 20 0'],...
        ['subregion 1  property 1 boundary 1 2 3 4 5 6']};
    clear adaptivity_options
    adaptivity_options.mesh_options =struct('quadratic',~true);
    adaptivity_options.targetnel=1650;
    adaptivity_options.convergence_rate=1.5;
    adaptivity_options.initial_mesh_size=5;
    adaptivity_options.nadapt=5;
    adaptivity_options.observer =@observer;
    
    function  model_data = make_model_data (fens,fes,groups,edge_fes,edge_groups)
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.conductivity =kappa;
        region.Q =Q;
        region.fes= fes;
        region.integration_rule =tri_rule(struct('npts',1));
        region.Rm =[];
        model_data.region{1} =region;
        
        clear essential
        essential.temperature=30;
        essential.node_list = [fenode_select(fens,struct('box',[48 48 0 48],...
            'inflate', 0.01)),fenode_select(fens,struct('box',[0 48 48 48],...
            'inflate', 0.01)),fenode_select(fens,struct('box',[0 20 0 13],...
            'inflate', 0.01))];
        model_data.boundary_conditions.essential{1} = essential;
    end
    
    model_data =heat_diffusion_adaptive_2D_steady_state(region_definition,@make_model_data,adaptivity_options);
    
    % Plot
    function  observer(step,model_data)
        disp(['Adaptive step ' num2str(step), ', ' num2str(count(model_data.region{1}.fes)) '  elements']);
        clear options
        model_data.postprocessing.z_scale = 1;
        model_data.postprocessing.camera =[-214.1067  -64.2844  201.6982   19.8741   23.1974   13.4603    0.5637    0.2108    0.7986 10.3396];
        model_data.postprocessing.draw_mesh = true;
        model_data =heat_diffusion_plot_raised_surface(model_data);
        pause(1);
    end
    %     if (graphics )
    %         if (Plot_errors)
    %             gv=graphic_viewer;
    %             gv=reset (gv,struct('limits', [0,48,0,48]));
    %             cmap=gray;cmap=cmap(end:-1:1,:);colormap(cmap);
    %             nvals=elerrs;
    %             dcm=data_colormap(struct ('range',[min(nvals),max(nvals)], 'colormap',cmap));%[min(nvals),max(nvals)]
    %             for i=1:count (gcells)
    %                 color = map_data (dcm, nvals(i));
    %                 draw(subset(gcells,i), gv, struct ('x',geom,'u',0*geom, 'facecolor',color));
    %             end
    %             view (2)
    %             lighting  none;
    %             nvalsrange=[min(nvals),max(nvals)];
    %             cbh=colorbar;
    %             set(cbh,...
    %                 'Position',[0.8 0.15 0.05 0.7],...
    %                 'YLim',[0,1],...
    %                 'YTick',[0,1],...
    %                 'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
    %             set(get(cbh,'XLabel'),'String',['Error']);
    %             pause(1);
    %         else
    %             gv=graphic_viewer;
    %             gv=reset (gv,[]);
    %             T=get(theta,'values');
    %             dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',jet));
    %             colorfield=field(struct ('name', ['colorfield'], 'data',...
    %                 map_data(dcm, T)));
    %             geomT=field(struct ('name', ['geomT'], ...
    %                 'data',[get(geom,'values'), get(theta,'values')]));
    %             draw(gcells, gv, struct ('x',geomT, 'u',0*geomT,...
    %                 'colorfield',colorfield, 'shrink',1));
    %             draw(gcells, gv, struct ('x',geom, 'u',0*geom, ...
    %                 'facecolor','none'));
    %             % draw_integration_points(femm, gv, struct ('x',geom,'u',0*geom, ...
    %             %     'theta', theta, 'scale', 35,'color','red'));
    %             xlabel('X [m]')
    %             ylabel('Y [m]')
    %             zlabel('Temperature [degrees C]')
    %             headlight(gv);
    %             camset (gv, [-315.9949  -97.6724  226.3628   17.7580   27.1126   12.2654         0         0    1.0000    7.6699])
    %         end
    %     end
    %
    %     [fens,gcells,groups,edge_gcells,edge_groups] ...
    %         = targe2_mesher_adapt(Region_definition,conn,x,hests,1.0,Mesh_options);
    
end