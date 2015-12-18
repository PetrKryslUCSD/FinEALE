function  THICK_ss_plate_pl_had
    disp('Simply-supported circular plate with pressure load');
    E=10.92e6;% The Young's modulus is taken  much larger to limit the amount of deflection
    nu=0.3;
    Magnitude=1.0;
    R=5.0;
    thickness=  1.0;
    Pressure=1.0;% Applied pressure
    
    Region_definition={...
        ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(R) ' ' num2str(0) ],...
        ['curve 2 line ' num2str(R) ' ' num2str(0) ' ' num2str(R) ' ' num2str(thickness) ],...
        ['curve 3 line ' num2str(R) ' ' num2str(thickness)  ' ' num2str(0) ' ' num2str(thickness) ],...
        ['curve 4 line ' num2str(0) ' ' num2str(thickness)  ' ' num2str(0) ' ' num2str(0) ],...
        ['subregion 1  property 1 boundary 1 2 3 4 ']};
    Mesh_options =struct('axisymm',true,'quadratic',true);
    Targetnel=600;
    convergence_rate=1.5;
    graphics=~false;
    sigj=4;
    
    % Mesh'
    mesh_size=R/4;
    [fens,fes,groups,edge_fes,edge_groups]...
        =targe2_mesher(cat(2,Region_definition,{['m-ctl-point constant ' num2str(mesh_size)]}), 1.0, Mesh_options);
    %     drawmesh({fens,fes},'fes','facecolor','red');  view(2); return
    
    
    for Adapt=1:7
        
        count(fes)
        
        % Material
        prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
        mater = material_deformation_linear_biax (struct('property',prop, ...
            'reduction','axisymm'));
        % Finite element block
        femm = femm_deformation_linear(struct ('material',mater, 'fes',fes,...
            'integration_rule',tri_rule (struct('npts',6))));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
        % Define the displacement field
        u   = clone(geom,'u');
        u   = u*0; % zero out
        % Apply EBC's
        ebc_fenids=fenode_select (fens,struct('box',[R R 0 thickness],'inflate',thickness/10000));
        ebc_fixed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        % Number equations
        u   = numberdofs (u);
        % Assemble the system matrix
        K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
        % Load
        efemm = femm_deformation_linear (struct ('material',mater, ...
            'fes',subset(edge_fes,edge_groups{3}),...
            'integration_rule',gauss_rule (struct('dim',1,'order', 2))));
        fi=force_intensity(struct('magn',[0;Pressure]));
        F = distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        % Solve
        u = scatter_sysvec(u, K\F);
        % get(u,'values')
        
        
        nodal_stress = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', 1:4);
        elerrs = flux_L2_error (femm, geom, u, [], nodal_stress);
        total_err=sqrt(sum(elerrs.^2));
        targeterr=sqrt(total_err^2/Targetnel);
        
        [hcurs, hests] =T3_mesh_sizes(fes.conn,geom.values,targeterr,elerrs,convergence_rate);
        
        fld1 = field_from_integration_points(femm, geom, u, [], 'Cauchy',sigj);
        nvals=fld1.values;%min(nvals),max(nvals)
        
        if (graphics)
            gv=graphic_viewer;
            gv=reset (gv,struct ('limits',[0 1.06*R -0.5*thickness 1.6*thickness]));
            scale=1e4;
            cmap = jet;
            nvalsrange=[min(nvals),max(nvals)];
            dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
            geom3=nodal_field(struct ('name', ['geom3'], ...
                'data',[geom.values,0.*nvals]));
            u3=nodal_field(struct ('name', ['u3'], ...
                'data',[u.values,0*nvals]));
            draw(femm,gv, struct ('x', geom3, 'u', +scale*u3,'colorfield',colorfield, 'edgecolor', 'black','shrink',1.0));
            %     draw(femm,gv, struct ('x', geom3, 'u', +0*u3,'facecolor','none', 'edgecolor','black'));
            draw_colorbar(gv, struct('colormap',cmap,'position',[0.85 0.15 0.05 0.7],...
                'minmax',nvalsrange,'label',['\sigma_{' num2str(sigj) '}']));
            pause(1);
        else
            plot(nvals(ebc_fenids),geom.values(ebc_fenids,2),'+r','Linewidth',3)
            if (sigj==1) xlabel '\sigma_{r}'
            elseif (sigj==2) xlabel '\sigma_{z}'
            elseif (sigj==3) xlabel '\sigma_{\theta}'
            elseif (sigj==4) xlabel '\sigma_{rz}'
            end
            ylabel 'z'
            pause(0.1);
        end
        
        
        [fens,fes,groups,edge_fes,edge_groups] ...
            = targe2_mesher_adapt(Region_definition,fes.conn,geom.values,hests,1.0,Mesh_options);
    end
    
    function idat=idi(idat, out, xyz, u, pc)
        [V,D]=eig(stress_4v_to_3x3t (mater,out([1, 2, 4, 3])));
        sD=sort(diag(D),'descend');
        if (sD(1)>=maxP1)
            maxP1=sD(1);
            maxP2=sD(2);
            maxP3=sD(3);
            maxxy=xyz;
            stress_4v_to_3x3t (mater,out);
        end
    end
    function idat=idi2d(idat, out, xyz, u, pc)
        [V,D]=eig(stress_3v_to_2x2t (mater,out([1, 2, 4])));
        sD=sort(diag(D),'descend');
        if (sD(1)>=maxP12d)
            maxP12d=sD(1);
            maxP22d=sD(2);
            maxxy2d=xyz;
            %              stress_3v_to_2x2t (mater,out([1, 2, 4]))
        end
    end
    
end%
