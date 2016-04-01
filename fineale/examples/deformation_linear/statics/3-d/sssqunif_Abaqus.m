function sssqunif
    
    % Simply-supported square plate with uniform distributed load.
    % MSC.Marc: simply supported square plate with uniform load
    % Demonstration problem 2.13.
    % Analytical solution from Roark's formulas, seventh edition, Table 11.4 (page 502).
    disp('% MSC.Marc: simply supported square plate with uniform load');
    % Parameters:
    um=physical_units_machine;;
    q=1*um('PSI'); a=60*um('in'); E=2e7*um('PSI'); h=00.3*um('in'); nu=0.3;
    w_max = 0.0444*q*a^4/(E*h^3);
    scale = a/4/w_max;
    graphics = false; % graphic output
    export_to_abaqus= ~true; SurfElType ='SFM3D4';;
    ElType ='C3D8';
    %     ElType ='C3D8H';
    %     ElType ='C3D8I';
    %     ElType ='C3D8IH';
    %          ElType ='C3D8RH';
    %         ElType ='C3D20R'; SurfElType ='SFM3D8';;
        % ElType ='C3D20H'; SurfElType ='SFM3D8';;


    % Mesh
    neqs=[]; normalized_deflections = [];
    for nt=2
        for na =2:2:8
            % tic;
            [fens,fes] = H8_block (a/2,a/2,h,na,na,nt);
            ir=gauss_rule(struct('dim',3,'order',2)); description  ='H8';
                                    [fens,fes] = H8_to_H20(fens,fes);
                                    ir=gauss_rule(struct('dim',3,'order',2)); description  ='H20R';
            %                         ir=gauss_rule(struct('dim',3,'order',3)); description  ='H20';
            % Material
            prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
            mater = material_deformation_linear_triax (struct('property',prop));
            % Finite element block
            femm = femm_deformation_linear(struct('material',mater,'fes',fes,...
                'integration_rule',ir));
            %             femm = femm_deformation_nonlinear_h8msgso(struct('material',mater,'fes',fes,...
            %                 'integration_rule',ir));
            % Geometry
            geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
            % Define the displacement field
            u   = 0*geom; % zero out
            % Specify the EBC's
            % Plane of symmetry
            ebc_fenids=fenode_select (fens,struct('box',[0 0 0 a/2 0 h],'inflate',h/100));
            u   = set_ebc(u, ebc_fenids, true, 1, 0.0);
            ebc_fenids=fenode_select (fens,struct('box',[0 a/2 0 0 0 h],'inflate',h/100));
            u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
            % Simple support
            ebc_fenids=fenode_select (fens,struct('box',[0 a/2 a/2 a/2 0 h],'inflate',h/100));
            u   = set_ebc(u, ebc_fenids, true, 3, 0.0);
            ebc_fenids=fenode_select (fens,struct('box',[a/2 a/2 0 a/2 0 h],'inflate',h/100));
            u   = set_ebc(u, ebc_fenids, true, 3, 0.0);
            % Apply
            u   = apply_ebc (u);
            % Number equations
            u   = numberdofs (u);
            % Assemble the system matrix
            %             femm  =update(femm,geom,u,u);
            K = stiffness(femm, sysmat_assembler_sparse,   geom, u);
            % Load
            bfe= mesh_boundary (fes,[]);
            ttopfl =fe_select (fens,bfe,struct('box',[-inf inf -inf inf h h],'inflate',h/100));
            ttopf=subset(bfe,ttopfl);
            fi=force_intensity(struct('magn',[0;0;-q]));
            qfemm = femm_deformation_linear(struct('material',mater,'fes',ttopf,...
                'integration_rule',gauss_rule(struct('dim',2,'order',3))));
            F = distrib_loads(qfemm, sysvec_assembler, geom, u, fi, 2);
            
            % Solve
            u = scatter_sysvec(u, K\F);
            corn=fenode_select (fens,struct('box',[0 0 0 0 0 h],'inflate',h/100));
            ucorn= gather_values(u, corn);
            nd=abs(mean(ucorn(:,3))/w_max);
            
            if (export_to_abaqus)
                date_now = clock;
                s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
                AE=Abaqus_exporter;
                AE.open([mfilename '-' ElType '-' num2str(na) '-' s '.inp']);
                AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'na=' num2str(na)]);
                AE.PART('PART1');
                AE.END_PART();
                AE.ASSEMBLY('ASSEM1');
                AE.INSTANCE('INSTNC1','PART1');
                AE.NODE(geom.values);
                AE.ELEMENT(ElType,'All',1,fes.conn);
                AE.ELEMENT(SurfElType,'Traction',count(fes)+1,ttopf.conn);
                AE.NSET_NSET('Corner',corn);
                AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
                AE.SOLID_SECTION('Material','Global','All','Hourglass');
                AE.SURFACE_SECTION('Traction');
                AE.NSET_NSET('xfix',find(u.is_fixed(:,1)));
                AE.NSET_NSET('yfix',find(u.is_fixed(:,2)));
                AE.NSET_NSET('zfix',find(u.is_fixed(:,3)));
                AE.END_INSTANCE();
                AE.END_ASSEMBLY();
                AE.MATERIAL('Material');
                AE.ELASTIC_ISOTROPIC(E,nu);
                %                                 AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
                AE.STEP_PERTURBATION('Linear solution');
                AE.DLOAD('ASSEM1.INSTNC1.Traction',[0;0;-q]);
                AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
                AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
                AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
                AE.NODE_PRINT('ASSEM1.INSTNC1.Corner');
                AE.ENERGY_PRINT();
                AE.END_STEP();
                AE.close();
                abaqus_job([AE.filename ]);
                AW=Abaqus_lck_watcher();
                AW.wait([AE.filename '.lck']);
                try
                    d= extract_displacement_from_abaqus_dat([AE.filename '.dat'],...
                        'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_CORNER',...
                        length(corn));
                catch,
                    d=[0,0, 0]; energy = 0;
                end
                nd=abs(mean(d(:,3))/w_max);
            
            end
            
            disp(['nel=',num2str(na),', Normalized deflection =', num2str(nd),' (deflection =', num2str(w_max/um('in')),')' ])
            %toc
            normalized_deflections = [normalized_deflections,nd];
            neqs = [neqs u.nfreedofs];
            if graphics
                gv=graphic_viewer;
                gv=reset (gv,[]);
                draw(fes, gv, struct ('x',1/um('in')*geom,'u', scale/um('in')*u, 'facecolor','none'));
                %             camset(gv, [-3.3445,-4.2190,3.1679,0.3162,0.5517, -0.3039 , 0, 0,1.0000,10.3396]);
                labels
            end
        end
    end
    if ~graphics
        % plot(neqs, normalized_deflections,'gd-','linewidth',3); hold on
        % plot(neqs, normalized_deflections,'ro-','linewidth',3); hold on
        % plot(neqs, normalized_deflections,'bx-','linewidth',3); hold on
        plot(neqs, normalized_deflections,name_to_style(description ),'linewidth',3); hold on
        xlabel ('Number of dofs')
        % ylabel ('Est. True Error of Deflection')
        ylabel ('Normalized deflection at center')
        %         % legend(Description,'Location','Southeast')
                set(gca,'linewidth',2);
        %         set(gca,'Position',[0.2 0.17 0.75 0.78]);
        options.FontSize=16;
        set_pub_defaults(gcf,options);
        grid on
        %         title([ 'Span to thickness ratio =' num2str(a/h) ])
    end