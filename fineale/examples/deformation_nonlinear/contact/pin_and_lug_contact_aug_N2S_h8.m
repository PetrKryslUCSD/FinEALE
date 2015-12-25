function  spherical_punch_aug_h20
    disp('Lug and pin: the cylindrical pin surface is  rigid');

    c =clock; disp(['Started ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
     % Parameters:
    U=physical_units_struct;
    E=210000*U.MEGA*U.PA;
    nu=0.3;
    L= 0.03*U.M; % in-plane dimension
    W = 0.03*U.M; % in-plane dimension
    a= 0.015*U.M; % hole radius
    H = 0.005*U.M; % thickness of the plate
    nL=5;nH=4;nW=5;na=15;tolerance=min([L,W,a,H])/1e5;
    tol = a*10e-7;
    sigma0=1*U.MEGA*U.PA;
    sscale=2000.0;
    penalty = 1*E*a*H;
    utol=1e-9*U.M;
    augtol=1e-3;
    Rrat= 0.9;
    surface_data.center=[(Rrat-1.0)*a,0,0];
    surface_data.R=Rrat*a;

    function  [penetration,normal] = get_penetration (surface_data,X,U)
        normal = get_normal(surface_data,X);
        p=X-surface_data.center;
        penetration = surface_data.R - dot(p,normal) - dot(U,normal);
    end
    function normal = get_normal(surface_data,x)
        x(3)=0;
        normal=(x-surface_data.center);
        normal= normal/norm(normal);
    end

    % Mesh
    [fens,fes]=Q4_elliphole(a,a,L,W,nL,nW,na,[]);
    [fens,fes] = H8_extrude_Q4(fens,fes,nH,@(x,i)([x,0]+[0,0,H/nH*i]));
    [fens1,fes1] = mirror_mesh(fens, fes,...
        [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
    [fens,fes1,fes2] = merge_meshes(fens, fes, fens1,    fes1, tolerance);
    fes=cat(fes1,fes2);
    
    bdry_fes = mesh_boundary(fes,  []);
    pullfacelist =fe_select(fens,bdry_fes,struct ('box',[L L -inf inf -inf inf],'inflate',tolerance));
    holefacelist =fe_select(fens,bdry_fes,struct ('cylinder',[0 0 0 a 20*H],'inflate',tolerance));
    %     clampfacelist =fe_select(fens,bdry_fes,struct ('box',[-L -L -inf inf -inf inf],'inflate',tolerance));
    cn = connected_nodes(subset(bdry_fes,holefacelist));
    penetration_fes=fe_set_P1 (struct('conn',cn));
    cn = connected_nodes(subset(bdry_fes,pullfacelist));
    Spring_fes=fe_set_P1 (struct('conn',cn));
    %     gv=drawmesh({fens,bdry_fes},'fes', 'facecolor','r')
    % gv=drawmesh({fens,subset(bdry_fes,holefacelist)},'gv',gv,'fes', 'facecolor','y')



    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    % Finite element block
    femm = femm_deformation_linear_h8msgso(struct('fes',fes, 'material',mater,...
            'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    penetrationfemm = femm_deformation_linear_penetration (struct('fes',penetration_fes,...
         'integration_rule', point_rule,'get_penetration',@get_penetration,'surface_data',surface_data,'penalty',penalty));
    Springfemm = femm_deformation_spring_grounded (struct('fes',Spring_fes,...
        'translation_stiffness_matrix',penalty/100000*diag([1,0,0])));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = clone(geom,'u');
    u   = u*0; % zero out
    % Apply EBC's
    ebc_fenids=fenode_select (fens,struct('box',[-inf inf -inf inf 0 0],'inflate',tolerance));
    u   = set_ebc(u, ebc_fenids, true, 3, 0.0);
    ebc_fenids=fenode_select (fens,struct('box',[-inf inf 0 0 -inf inf],'inflate',tolerance));
    u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
    %     ebc_fenids=fenode_select (fens,struct('box',[-a -a 0 0 0 0],'inflate',tolerance));
    %     u   = set_ebc(u, ebc_fenids, true, 1, 0.0);
    %     ebc_fenids=connected_nodes(subset(bdry_fes,clampfacelist));
        %     u   = set_ebc(u, ebc_fenids, true, [], 0.0);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
        
    
    u   = u*0; % zero out
    
    tic;
    Ks = stiffness(femm, sysmat_assembler_sparse, geom, u);
    Ksp =stiffness(Springfemm, sysmat_assembler_sparse, geom, u);
    disp(['Stiffness matrix assembly ' num2str(toc)])
    % Load
    efemm = femm_deformation_linear (struct ('material',mater, ...
        'fes',subset(bdry_fes,pullfacelist),...
        'integration_rule',gauss_rule (struct('dim',2,'order', 2))));
    fi=force_intensity(struct('magn',[sigma0;0;0]));
    Fl = distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
    R = resultant_force(u,Fl)
    %     us = K\(F);
    %     u = scatter_sysvec(u, us);
    %     gl = penetration_loads(penetrationfemm, sysvec_assembler, geom, u);
    %     resultant_force(u,gl)
    lm =0*Fl;
    du=0*u;
    u1=0*u;
            gv=reset (graphic_viewer,[]);
    for augit= 1:20
        disp([' Augmented Iteration ' num2str(augit)])
        iter=1;
        while 1
            gl = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
            %             resultant_force(u,gl)
            U1= gather_sysvec(u1);
            Fr  = -(Ksp*U1+Ks*U1);
            Kg =stiffness(penetrationfemm, sysmat_assembler_sparse, geom, u1);
            K = Ks + Ksp + Kg;
            du = scatter_sysvec(du, K\(-lm+gl+Fr+Fl));
            u1 = u1 + du;   % increment displacement
            disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
            if (max(abs(du.values)) < utol) break; end;                    % convergence check
            iter=iter+1;
        end
        lm=lm-gl;
        disp(['   Augmented Lagrange It. ' num2str(augit) ': ||gl||/||lm||=' num2str(norm(gl)/norm(lm))]);
        if (norm(gl)<=augtol*norm(lm)), break,end
                 draw(bdry_fes,gv, struct ('x', geom, 'u', sscale*u1,'facecolor','red'));
             draw(subset(bdry_fes,holefacelist),gv, struct ('x', geom, 'u', (0)*u1,'facecolor','y'));
    end
    u=u1;% Now we have our final solution
    %         if isempty(progressbar((step/nsteps),h)), return, end;
    
    c =clock; disp(['Stopped ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
    % Plot
    gv=graphic_viewer;
    gv=reset (gv,[]);
    scale=1;
    %         draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','none'));
    draw(bdry_fes,gv, struct ('x', geom, 'u', sscale*u,'facecolor','red'));
     draw(subset(bdry_fes,holefacelist),gv, struct ('x', geom, 'u', (0)*u,'facecolor','y'));
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on
end
