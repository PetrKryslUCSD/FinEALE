function squash_thecube(ntmult)
    if ~exist('ntmult','var'), ntmult=1; end
    % A new locking-free brick element technique for large
    % deformation problems in elasticity p
    % S. Reese a, *, P. Wriggers a , B.D. Reddy
    %     kappa=400889.8;
    mu=80.194;
    %     %     [x,y]=solve ('400889.8-x/3/( 1-2*y)','80.194-x/2/(1+y)')
    %     E = 240.56595906096809707175008200418;
    %     nu=0.49989998666339188138607677634352;
    lambda =400889.8;;%lambda = E * nu / (1 + nu) / (1 - 2*(nu));
    %     [x,y]=solve ('400889.8-x*y/( 1+y)/(1-2*y)','80.194-x/2/(1+y)')
    E = 240.56596;
    nu=0.4999;;
    
    %     nu=.0;
    Fmag= 80*4;
    utol = 1e-9;
    nincr = 40;
    graphics = ~false;
    scale=1;
    nt = [ntmult*[2,2,2]];
    %     nt = [4*[2,2,2]];
    
    prop = property_deformation_neohookean (struct('E',E,'nu',nu));
    mater = material_deformation_neohookean_triax(struct('property',prop));
    alp=min([1.0,max([(1-2*nu)*5,0.00001])]);;
    nuhat =1/2-(1-2*nu)/2/alp;%A reasonable choice.   No locking anyway.
    stabprop = property_deformation_neohookean (struct('E',E,'nu',nuhat));
    stabmater = material_deformation_neohookean_triax(struct('property',stabprop));
    
    %
    function [fens,fes] = h8T4(nt)
        [fens,fes]=t4blocka(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
    end
    function [fens,fes] = h8T4d(nt)
        [fens,fes]=t4blocka(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
        [fens1,fes1]=t4blockb(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, 1/1000);
        fes=cat(2,fes1,fes2);
    end
    function [fens,fes] = h8H8(nt)
        [fens,fes]=H8_block(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
    end
    function [fens,fes] = h8H8u(nt)
        [fens,fes]=H8_block_u(1/2, 1/2, 1/2, nt(1)/2,nt(2)/2,nt(3)/2);
        [fens1,fes1]=mirror_mesh(fens,fes,[-1,0,0],[0,0,0],@(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, 1/1000);
        fes=cat(fes1,fes2);
        [fens1,fes1]=mirror_mesh(fens,fes,[0,-1,0],[0,0,0],@(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, 1/1000);
        fes=cat(fes1,fes2);
        [fens1,fes1]=mirror_mesh(fens,fes,[0,0,-1],[0,0,0],@(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, 1/1000);
        fes=cat(fes1,fes2);
        fens=translate_mesh(fens, [0.5, 0.5, 0.5]);
        count(fes)
        disp(' ')
        %         gv=drawmesh( {fens,fes},'facecolor','r');
    end
    function [fens,fes] = h8H27(nt)
        [fens,fes]=H8_block(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
        [fens,fes] = H8_to_H27(fens,fes);
    end
    function [fens,fes] = h8H64(nt)
        [fens,fes]=block(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
        [fens,fes] = H8_to_H64(fens,fes);
    end
    
    eix=1;
    clear eltyd
    
    %     eltyd(eix).description='H8MSGSO(U)';
    %     eltyd(eix).mf =@h8H8u;
    %     eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgso(struct ('material',mater, ...
    %         'fes',fes,  'stabilization_material',stabmater,...
    %         'integration_rule',gauss_rule(struct('dim',3,'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
    %     eix=eix+1;
    
        eltyd(eix).description='H8MSGSO(U)';
        eltyd(eix).mf =@h8H8u;
        eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
            'integration_rule',gauss_rule(struct('dim',3,'order',2))));
        eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
        eix=eix+1;
    
    %     eltyd(eix).description='H8MSGS(U)';
    %     eltyd(eix).mf =@h8H8u;
    %     eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgs(struct ('material',mater, 'fes',fes, ...
    %         'integration_rule',gauss_rule(struct('dim',3,'order',2)),'stab_fraction',0.75));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
    %     eix=eix+1;
    
    %     eltyd(eix).description='H8MSGS';
    %     eltyd(eix).mf =@h8H8;
    %     eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgs(struct ('material',mater, 'fes',fes, ...
    %         'integration_rule',gauss_rule(struct('dim',3,'order',2)),'stab_fraction',0.75));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
    %     eix=eix+1;
    
    % eltyd(eix).label ='T4';
    % eltyd(eix).mf =@h8T4;
    % eltyd(eix).blf =@feblock_defor_ss;
    % eltyd(eix).integration_rule= tet_rule(1);
    % eltyd(eix).surface_integration_rule=tri_rule(1);
    % eltyd(eix).styl='rv-.';
    % eix=eix+1;
    %
    % eltyd(eix).label ='H8';
    % eltyd(eix).mf =@h8H8;
    % eltyd(eix).blf =@feblock_defor_nonlinear;
    % eltyd(eix).integration_rule= gauss_rule (3,2);
    % eltyd(eix).surface_integration_rule=gauss_rule (2,2);
    % eltyd(eix).styl='rx-.';
    % eix=eix+1;
    
    % eltyd(eix).label ='NICE-H8';
    % eltyd(eix).mf =@h8H8;
    % eltyd(eix).blf =@feblock_defor_nonlinear_nice;
    % eltyd(eix).integration_rule= tensprod_nq_rule (struct('dim',3,'order',1));
    % eltyd(eix).surface_integration_rule=tensprod_nq_rule(struct('dim',2,'order',1));
    % eltyd(eix).styl='rx-.';
    % eix=eix+1;
    
    %     eltyd(eix).label ='NICE-T4';
    %     eltyd(eix).mf =@h8T4;
    %     eltyd(eix).blf =@feblock_defor_nonlinear_nice;
    %     eltyd(eix).integration_rule= simplex_nq_rule (3);
    %     eltyd(eix).surface_integration_rule=tri_rule(1);
    %     eltyd(eix).styl='m^-';
    %     eix=eix+1;
    
    
    % eltyd(eix).label ='NICE-T4d';
    % eltyd(eix).mf =@h8T4d;
    % eltyd(eix).blf =@feblock_defor_nonlinear_nice;
    % eltyd(eix).integration_rule= simplex_nq_doubled_rule (3);
    % eltyd(eix).surface_integration_rule=simplex_nq_doubled_rule (2);
    % eltyd(eix).styl='m^-';
    % eix=eix+1;
    
    %     eltyd(eix).label ='H27';
    %     eltyd(eix).mf =@h8H27;
    %     eltyd(eix).blf =@femm_deformation_nonlinear;
    %     eltyd(eix).integration_rule=gauss_rule(struct('dim',3,'order',3));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',3));
    %     eltyd(eix).styl='ks-';
    %     eix=eix+1;
    
    % eltyd(eix).label ='NICE-H27';
    % eltyd(eix).mf =@h8H27;
    % eltyd(eix).blf =@feblock_defor_nonlinear_nice;
    % eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
    % eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
    % eltyd(eix).styl='mh-';
    % eix=eix+1;
    
    % eltyd(eix).label ='H20R';
    % eltyd(eix).mf =@uh8cylH20;
    % eltyd(eix).blf =@feblock_defor_ss;
    % eltyd(eix).integration_rule=gauss_rule (3,2);
    % eltyd(eix).surface_integration_rule=gauss_rule(2, 2);
    % eltyd(eix).styl='ro-';
    % eix=eix+1;
    % %
    
    % eltyd(eix).label ='NICE-H64';
    % eltyd(eix).mf =@h8H64;
    % eltyd(eix).blf =@feblock_defor_nonlinear_nice;
    % eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',3));
    % eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',3));
    % eltyd(eix).styl='mh-';
    % eix=eix+1;
    
    
    function [neqns,energy] =Compute(nt)
        %          Mesh
        [fens,fes] = eltyd(eix).mf (nt);
        %         [fens, fes]=renum_mesh(fens, fesAnd, 'symrcm');
        femm = eltyd(eix).blf (fes);
        %         [fens, fes]=renum_mesh_renumber(fens, fes,...
        %             renum_mesh_numbering(fens,genconn(feb), 'symrcm'));
        %         [fens, fes]=renum_mesh_renumber(renum_mesh_numbering(fens,genconn(feb), 'symamd'),fens, fes);
        bfes =mesh_boundary(fes, []);
        %         icl = fe_select(fens, bfes, struct('box', [0,1,0,1, 1,  1],'inflate',1/100));
        loadcl = fe_select(fens, bfes, struct('box', [0,0.5,0,0.5, 1,  1],'inflate',1/100));
        cncl = fenode_select(fens, struct('box', [0,0,0,0,1,1],'inflate',1/100));
        %          gv=drawmesh( {fens,fes},'facecolor','none'); gv=drawmesh( {fens,bfes(icl)},'gv',gv,'facecolor','blue');
        %         view(3); pause(1)
        % Material
        % Finite element block
        efemm = femm_deformation (struct ('material',mater, 'fes',subset(bfes,loadcl),...
            'integration_rule',eltyd(eix).surface_integration_rule));
        sfemm = femm_deformation (struct ('material',mater, 'fes',bfes,...
            'integration_rule',eltyd(eix).surface_integration_rule));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = 0*geom; % zero out
        % Apply EBC's
        %         Bottom surface
        ebc_fenids=fenode_select (fens,struct('box',[0,1,0,1,0,0],'inflate',1/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+3;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        %         Top surface
        ebc_fenids=fenode_select (fens,struct('box',[0,1,0,1,1,1],'inflate',1/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+1;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        % Symmetry planes
        ebc_fenids=fenode_select (fens,struct('box',[0,0,0,1,0,1],'inflate',1/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+1;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[0,1,0,0,0,1],'inflate',1/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        
        % Number equations
        u   = numberdofs (u);
        % Now comes the nonlinear solution
        tup = 1;
        u = u*0; % zero out the displacement
        utol =         utol*u.nfreedofs;
        us={};
        
        if (graphics),
            gv=reset(clear(graphic_viewer,[]),[]);
            cmap = jet;
        end
        
        t=0; % time=load magnitude
        femm  =associate_geometry(femm,geom);
        incr=1; dt=tup / nincr;
        while (incr <= nincr)
            t = t + dt;
            disp(['Increment ' num2str(incr) ]); % pause
            % Initialization
            u1 = u; % guess
            u1 = apply_ebc(u1);
            du = 0*u; % this will hold displacement increment
            du = apply_ebc(du);
            iter=1;
            %             gv=reset(clear(gv,[]),[]);
            %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
            femm1=femm;
            while 1
                Load=zeros(3,1); Load(3) =-Fmag*t;
                fi=force_intensity(struct('magn',Load));
                FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
                F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u,dt);       % Internal forces
                K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u,dt) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u,dt);
                % Displacement increment
                du = scatter_sysvec(du, K\F);
                R0 = dot(F,gather_sysvec(du));
                F = FL + restoring_force(femm1,sysvec_assembler, geom,u1+du,u,dt);       % Internal forces
                R1 = dot(F,gather_sysvec(du));
                a = R0/R1;
                if ( a<0 )
                    eta = a/2 +sqrt((a/2)^2 -a);
                else
                    eta =a/2;
                end
                if (imag(eta)~=0)
                    disp('######################  Inverted elements?')
                end
                eta=min( [eta, 1.0] );
                u1 = u1 + eta*du;   % increment displacement
                if graphics, pause(0.5); Cam =camget(gv); end
                disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
                if (max(abs(du.values)) < utol) break; end;                    % convergence check
                iter=iter+1;
            end
            [~,femm]  =restoring_force(femm1,sysvec_assembler, geom,u1,u,dt); 
            disp(['    Converged for t=' num2str(t)]); % pause
            u = u1;                                               % update the displacement
            if graphics
                gv=reset(clear(gv,[]),[]);
                camset (gv,Cam);
                %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
                %             draw(sfemm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
                %                 fld = field_from_integration_points(femm, geom, u1, [], 'pressure',1);
                %                 nvals=get(fld,'values');%min(nvals),max(nvals)
                %                 nvalsrange=[min(nvals),max(nvals)];
                %                 dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
                %                 colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
                %                 draw(sfemm,gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink',1.0));
                draw(femm,gv, struct ('x', geom, 'u', u,'facecolor','none', 'shrink',1.0));
                % draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
                %                 colormap(cmap);
                %                 cbh=colorbar;
                %                 set(cbh,...
                %                     'Position',[0.815 0.15 0.05 0.7],...
                %                     'YLim',[0,1],...
                %                     'YTick',[0,1],...
                %                     'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
                %                 set(get(cbh,'XLabel'),'String','pressure');
                pause(0.5); Cam =camget(gv);
            end
            us{end+1} =u;
            utip=gather_values(u,cncl);
            disp(['    Center compression =' num2str(-utip(3))]); % pause
            incr = incr + 1;
        end
        clear K
        save ([mfilename  '-nt' num2str(nt) '-' eltyd(eix).description])
    end
    
    for eix = 1:length(eltyd)
        Compute(nt);
    end
end