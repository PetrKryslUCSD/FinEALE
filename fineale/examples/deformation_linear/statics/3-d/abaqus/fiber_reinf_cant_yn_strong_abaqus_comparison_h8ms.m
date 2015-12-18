% Fiber-reinforced cantilever.
function fiber_reinf_cant_yn_strong_abaqus_comparison_h8ms
    u= physical_units_struct;
    % Anisotropic, but much stronger anisotropy
    E1=100000*1e9*u.PA; E2=1e9*u.PA; E3=E2; G12=0.2e9*u.PA;  G13=G12; G23=0.2e9*u.PA;
    nu12= 0.25; nu13= 0.25; nu23= 0.25;
    aangle =-45;
    % [xestim, beta, c, residual] = richextrapol(uz(end-2:end),hs)
    uz_ref =  -5.9189e-06;
    prop = property_deformation_linear_ortho (...
        struct('E1',E1,'E2',E2,'E3',E3,...
        'G12',G12,'G13',G13,'G23',G23,...
        'nu12',nu12,'nu13',nu13,'nu23',nu23));
    
    a=90*u.MM; b=10*u.MM;  t=20*u.MM;
    q0 = -1000*u.PA;
    
    [na,nb,nt] =adeal(1*[1,1,1]);
    refs = [2,3,4,6];
    %         refs = [1,2,3];
    %            refs = [2,4,8];
    tolerance =t/nt/100;
    u_scale  = 10000;
    graphics = ~true;
    export_to_abaqus=  true; SurfElType ='SFM3D4';;
    %     ElType ='C3D8IH';
    ElType ='C3D8';
        ElType ='C3D8H';
                ElType ='C3D8I';
                ElType ='C3D8RH';
        %         ElType ='C3D20R'; SurfElType ='SFM3D8';;
    
    Rm= rotmat(aangle/180*pi* [0,1,0]);
    
    mater = material_deformation_linear_triax (struct('property',prop ));
    
prop_stab=property_deformation_linear_iso(struct('E',00.2*G23,'nu',0));
mater_stab = material_deformation_linear_triax (struct('property',prop_stab ));

    clear eltyd
    eix=1;

        
    % eltyd(eix).description ='H8MSNIS';
    % eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    % eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msnis(struct('fes',fes,'material',mater,'Rm',Rm,'nconstrained',1, 'stab_fraction',0.15));
    % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    % eix=eix+1;
        
          eltyd(eix).description ='H8MSGSO';
    eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
        'Rm',Rm,'integration_rule',gauss_rule(struct('dim',3,'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eix=eix+1;
    

    %           eltyd(eix).description ='H8MSGS';
    %     eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,'material',mater,'Rm',Rm,'integration_rule',gauss_rule(struct('dim',3,'order',2)), 'stab_fraction',0.1,'weights',E2/E1));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eix=eix+1;
    %
    %
    %           eltyd(eix).description ='H8MSGS';
    %     eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,...
    %         'material',mater,'Rm',Rm,'integration_rule',gauss_rule(struct('dim',3,'order',2)), 'stabilization_material',mater_stab));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eix=eix+1;
    
    %         eltyd(eix).description ='H20R';
    %         eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'Rm',Rm,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %         eix=eix+1;
    
    
    %     eltyd(eix).description ='H20-Bbar';
    %     eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
    %         'Rm',Rm,...
    %        'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %         'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
    %         'pv_bfun',@(p)[1;p(1);p(2);p(3)],'nconstrained',1));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eix=eix+1;
    %
    %             eltyd(eix).description ='H8';% tetrahedron
    %             eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %             eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
    %                 'Rm',Rm,...
    %                 'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %             eix=eix+1;
       
    %     eltyd(eix).description ='H8-BbarX';
    %     eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbarx(struct('fes',fes,'material',mater,...
    %         'Rm',Rm,...
    %         'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %         'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    %         'pv_bfun',[],'nconstrained',6,'psi',1-1e-7));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eix=eix+1;
    
    %         eltyd(eix).description ='H8-Bbar';
    %         eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
    %             'Rm',Rm,...
    %             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    %             'pv_bfun',@(p)[polynomial_basis(3,0,p)],'nconstrained',1));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %         eix=eix+1;
    
    % eltyd(eix).description ='H8-BbarX';
    % eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    % eltyd(eix).femmf =@(fes)femm_deformation_linear_bbarx(struct('fes',fes,'material',mater,...
    %     'Rm',Rm,...
    %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    %     'pv_bfun',[],'nconstrained',3));
    % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    % eix=eix+1;
    
    % eltyd(eix).description ='H8-BbarX';
    % eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    % eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
    %     'Rm',Rm,...
    %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    %     'pv_bfun',@(p)[polynomial_basis(3,0,p),0*polynomial_basis(3,0,p),0*polynomial_basis(3,0,p);...
    %     0*polynomial_basis(3,0,p),polynomial_basis(3,0,p),0*polynomial_basis(3,0,p);...
    %     0*polynomial_basis(3,0,p),0*polynomial_basis(3,0,p),polynomial_basis(3,0,p)]));
    % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    % eix=eix+1;
    
    % eltyd(eix).description ='H20R';
    % eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    % eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %     'Rm',Rm,...
    %     'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    % eix=eix+1;
    % %
    %
    eltyd(eix).description ='H20-Bbar';
    eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
        'Rm',Rm,...
       'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
        'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
        'pv_bfun',@(p)[1;p(1);p(2);p(3)],'nconstrained',1));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eix=eix+1;
    %
    %                eltyd(eix).description ='H8-Bbar-ISO';
    %             eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
    %             eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
    %                 'Rm',Rm,...
    %                 'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    %                 'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    %                 'nconstrained',[]));
    %             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %             eix=eix+1;
    
    legends ={};
    for eix = 1:length(eltyd)
        ns=[];
        uzs=[];
        uys=[];
        
        %         for ref = [4,8,16]
        for ref = refs
            
            % Mesh
            [fens,fes] = eltyd(eix).mf(ref);
            
            % Compose the model data
            clear model_data
            model_data.fens =fens;
            
            clear region
            region.fes= fes;
            region.femm= eltyd(eix).femmf(fes);
            model_data.region{1} =region;
            
            clear essential
            essential.component= [1,2,3];
            essential.fixed_value= 0;
            essential.node_list = fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
                'inflate',tolerance));
            model_data.boundary_conditions.essential{1} = essential;
            
            clear traction
            bdry_fes = mesh_boundary(fes, []);
            bcl = fe_select(fens, bdry_fes, ...
                struct ('box',[a,a,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
            traction.fes =subset(bdry_fes,bcl);;
            traction.traction= [0;0; q0];
            traction.integration_rule =eltyd(eix).surface_integration_rule;
            model_data.boundary_conditions.traction{1} = traction;
            
            model_data.renumber=true;
            
            % Solve
            model_data =deformation_linear_statics(model_data);
            
            
            nc=[fenode_select(fens, struct('box', [a,a,b,b,0,0],...
                'inflate',tolerance))];
            ab_displacement=gather_values(model_data.u,nc);
            abz =mean(ab_displacement(:,3) );
            aby =mean(ab_displacement(:,2) );
            
            if (export_to_abaqus)
                date_now = clock;
                s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
                AE=Abaqus_exporter;
                AE.open([mfilename '-' ElType '-' num2str(ref) '-' s '.inp']);
                AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'ref=' num2str(ref)]);
                AE.PART('part1');
                AE.END_PART();
                AE.ASSEMBLY('ASSEM1');
                AE.INSTANCE('INSTNC1','PART1');
                AE.NODE(model_data.geom.values);
                AE.ELEMENT(ElType,'All',1,fes.conn);
                AE.ELEMENT(SurfElType,'Traction',count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
                AE.NSET_NSET('Corner',nc);
                AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
                AE.ORIENTATION('Rm', Rm(:,1), Rm(:,2));
                AE.SOLID_SECTION('Material','Rm','All');
                AE.SURFACE_SECTION('Traction');
                AE.NSET_NSET('xfix',find(model_data.u.is_fixed(:,1)));
                AE.NSET_NSET('yfix',find(model_data.u.is_fixed(:,2)));
                AE.NSET_NSET('zfix',find(model_data.u.is_fixed(:,3)));
                AE.END_INSTANCE();
                AE.END_ASSEMBLY();
                AE.MATERIAL('Material');
                AE.ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
                AE.STEP_PERTURBATION('Linear solution');
                AE.DLOAD('ASSEM1.INSTNC1.Traction',model_data.boundary_conditions.traction{1}.traction);
                AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
                AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
                AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
                AE.NODE_PRINT('ASSEM1.INSTNC1.Corner');
                AE.ENERGY_PRINT();
                AE.END_STEP();
                AE.close();
                %                 delete([AE.filename '.dat']);
                system(['abaqus job=' [AE.filename ]]);
                AW=Abaqus_lck_watcher();
                AW.wait([AE.filename '.lck']);
                try
                    AR=Abaqus_dat_reader();
                    d= AR.extract_displacement_from_abaqus_dat([AE.filename '.dat'],...
                        'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_CORNER');
                catch,
                    d=[0,0, 0]; energy = 0;
                end
                abz=d(3);
                aby=d(2);
                eltyd(eix).description=ElType;
            end
            
            ns=[ns,model_data.u.nfreedofs];
            uzs=[uzs,abz];
            uys=[uys,aby];
            
            if (graphics )
                figure
                model_data.postprocessing.u_scale= u_scale;
                model_data.postprocessing.stress_component=2;
                %                     options=deformation_plot_stress(model_data)
                model_data=deformation_plot_deformation(model_data);
                %             anchor =[0,5*b,0]; tip = anchor +3*b*Rm(:, 1)';
                %             draw_arrow (model_data.postprocessing.gv,...
                %                 anchor, tip-anchor, struct('facecolor', 'black'));
                %             draw_polyline(model_data.postprocessing.gv,...
                %                 [anchor;  [tip(1),anchor(2),anchor(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
                %             draw_polyline(model_data.postprocessing.gv,...
                %                 [anchor;  [anchor(1),tip(2),anchor(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
                %             draw_polyline(model_data.postprocessing.gv,...
                %                 [anchor;  [anchor(1),anchor(2),tip(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
            end
        end
        format long
        disp(['Data{end+1}=[']);
        disp(num2str(refs));
        disp(num2str(ns));
        disp(num2str(uys,'%15.12e '));
        disp(num2str(uzs,'%15.12e '));
        disp(['];' 'Description{end+1}=''' eltyd(eix).description '''; %']);
        
        %         [xestim, beta, c, residual] = richextrapol(uzs(end-2:end),refs(end-2:end));
        if (~graphics )
            loglog(ns,abs(uzs/uz_ref-1),name_to_style(eltyd(eix).description),'linewidth',3); hold on
            %                 plot(ns,uys,eltyd(eix).styl,'linewidth',3); hold on
            %                 plot(diff(ns),diff(uzs),eltyd(eix).styl,'linewidth',3); hold on
            figure (gcf); pause (1)
        end
        legends{end+1} =eltyd(eix).description;
    end
    legend (legends);
    title (['Fiber dir =[' num2str(Rm(:,1)') ']'])
    labels(  'Number of equations', 'Estimated true error')
    grid on
    
    set_decades_on_axis (gca)
    set_pub_defaults
    % printeps(gcf,'motiv-fiber_reinf_cant_z_strong','fontsize',20)
    % clear options
    % options.u_scale= scale;
    % options.stress_component= 1;
    % options=deformation_plot_stress(model_data, options)
end
