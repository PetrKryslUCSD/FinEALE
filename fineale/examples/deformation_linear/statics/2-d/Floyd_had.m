% Floyd's pressure vessel example (from the book by Babuska, Szabo). H-adaptive analysis.
function  Floyd_had
    E= 1435;% PSI
    nu=0.49;
    % geometry
    R1= 6;%  inches
    Rf=0.15;% Inches
    R2=R1-0.6;% inches
    L1=3.15;%Inches
    L2=1.5;%Inches
    % fixed Traction
    Pressure=2.61;%PSI
    Region_definition={...
        ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(R1) ' ' num2str(0) ],...
        ['curve 2 line ' num2str(R1) ' ' num2str(0) ' ' num2str(R1) ' ' num2str(L1) ],...
        ['curve 3 line ' num2str(R1) ' ' num2str(L1) ' ' num2str(R2) ' ' num2str(L1) ],...
        ['curve 4 line ' num2str(R2) ' ' num2str(L1) ' ' num2str(R2) ' ' num2str(L2+Rf) ],...
        ['curve 5 arc ' num2str(R2) ' ' num2str(L2+Rf) ' ' num2str(R2-Rf) ' ' num2str(L2) ' center ' num2str(R2-Rf) ' ' num2str(L2+Rf) ],...
        ['curve 6 line ' num2str(R2-Rf) ' ' num2str(L2) ' ' num2str(0) ' ' num2str(L2) ],...
        ['curve 7 line ' num2str(0) ' ' num2str(L2) ' ' num2str(0) ' ' num2str(0) ],...
        ['subregion 1  property 1 boundary 1 2 3 4 5 6 7']};
    Mesh_options =struct('axisymm',true,'quadratic',true);
    Targetnel=600;
    convergence_rate=1.5;
    graphics=1;
    sigj=1;

    % Mesh'
    mesh_size=8/(2^3);
    [fens,fes,groups,edge_fes,edge_groups]...
        =targe2_mesher(cat(2,Region_definition,{['m-ctl-point constant ' num2str(mesh_size)]}), 1.0, Mesh_options);
%     drawmesh({fens,fes},'fes','facecolor','red');  view(2); return
    
    
    for Adapt=1:3
        
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
        ebc_fenids=fenode_select (fens,struct('box',[0 R1 L1 L1],'inflate',L1/1000));
        ebc_fixed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        % The axis of symmetry
        ebc_fenids=fenode_select (fens,struct('box',[0 0 0 L1+L2],'inflate',L1/1000));
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
        F = nz_ebc_loads(femm, sysvec_assembler, geom, u);
        efemm = femm_deformation_linear (struct ('material',mater, ...
            'fes',subset(edge_fes,edge_groups{4}),...
            'integration_rule',gauss_rule (struct('dim',1,'order', 2))));
        fi=force_intensity(struct('magn',[Pressure;0]));
        F = F + distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        efemm.fes = subset(edge_fes,edge_groups{5});
        fi=force_intensity(struct('magn',@(x) (Pressure*(x-[R2-Rf,L2+Rf])'/norm((x-[R2-Rf,L2+Rf])))));
        F = F + distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        efemm.fes = subset(edge_fes,edge_groups{6});
        fi=force_intensity(struct('magn',[0;-Pressure]));
        F = F + distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        % Solve
        u = scatter_sysvec(u, K\F);
        % get(u,'values')
        
        maxP1=-Inf; maxP2=-Inf; maxP3=-Inf;  maxxy=0;
        conns = fes.conn;
        for q=1:count(fes)
            conn = conns(q,:);
            x=gather_values(geom,conn);
            Cauchy = inspect_integration_points(femm, geom, u, [], ...
                q, struct ('output','Cauchy'), @idi,  []);
        end
        maxP12d=-Inf; maxP22d=-Inf;  maxxy2d=0;
        for q=1:count(fes)
            conn = conns(q,:);
            x=gather_values(geom,conn);
            Cauchy = inspect_integration_points(femm, geom, u, [], ...
                q, struct ('output','Cauchy'), @idi2d,  []);
        end
        disp([' Mesh size: ' num2str(mesh_size)]);
        disp([' Location 3d:  [' num2str(maxxy) '],  Principal stresses 3d: [' num2str(maxP1) ',' num2str(maxP2) ',' num2str(maxP3) ']' ]);
        disp([' Location 2d:  [' num2str(maxxy2d) '],  Principal stresses 2d: [' num2str(maxP12d) ',' num2str(maxP22d) ']' ]);
        %     return
        
        nodal_stress = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', 1:4);
        elerrs = flux_L2_error (femm, geom, u, [], nodal_stress);
        total_err=sqrt(sum(elerrs.^2));
        targeterr=sqrt(total_err^2/Targetnel);
        
        [hcurs, hests] =T3_mesh_sizes(fes.conn,geom.values,targeterr,elerrs,convergence_rate);
        
        if (graphics)
            gv=graphic_viewer;
            gv=reset (gv,struct ('limits',[0 1.06*R1 -0.5*L1 1.1*L1]));
            scale=5;
            cmap = jet;
            fld1 = field_from_integration_points(femm, geom, u, [], 'Cauchy',1);
            fld2 = field_from_integration_points(femm, geom, u, [], 'Cauchy',2);
            fld3 = field_from_integration_points(femm, geom, u, [], 'Cauchy',3);
            fld4 = field_from_integration_points(femm, geom, u, [], 'Cauchy',4);
            stresses = [fld1.values,fld2.values,fld3.values,fld4.values];
            Ps=zeros(size(stresses,1),3);
            maxP1=-Inf;
            for q= 1:size(stresses,1)
                [V,D]=eig(stress_3v_to_2x2t (mater,stresses(q,[1, 2, 4])));
                sD=sort(diag(D),'descend');
                Ps(q,1:2)=sD;Ps(q,3)=stresses(q,3);
            end
            nvals=Ps(:,sigj);%min(nvals),max(nvals)
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
    assignin('caller','fineale_test_passed',((norm(1.0e+002 *[1.026586000267639   0.104118128154772]-[maxP12d,maxP22d]))<1e-9));
end%
