% Fiber-reinforced cantilever plate.
function plate_2dir_strong_abaqus_comparison
    u= physical_units_struct;
    % Anisotropic, but much stronger anisotropy
    E1=1000*1e9*u.PA; E2=1000*1e9*u.PA; E3=1*1e9*u.PA; G12=0.2e9*u.PA;  G13=G12; G23=0.2e9*u.PA;
    nu12= 0.25; nu13= 0.25; nu23= 0.25;
    aangle =+45;
    % [xestim, beta, c, residual] = richextrapol([],[1,2,4])
    uz_ref = -1.779724511204966e-06; ;  % -2.144422740010290e-06;%-2.111838493841995e-06-2.258294496806001e-06
    prop = property_deformation_linear_ortho (...
        struct('E1',E1,'E2',E2,'E3',E3,...
        'G12',G12,'G13',G13,'G23',G23,...
        'nu12',nu12,'nu13',nu13,'nu23',nu23));
    
    a=90*u.MM; b=100*u.MM;  t=20*u.MM;
    q0 = -1000*u.PA;
    
    [na,nb,nt] =adeal(1*[1,1,1]);
    refs = [2,3,4,5,6,7,8];
        refs = [2,3,4,5];
    %         refs = [3,6,12];
    %             refs = [4,8,16,32]/2;
    %         refs = [8];
    %                     refs = [13];
    %                                        refs = [3];
    tolerance =t/nt/100;
    u_scale  = 5000;
    graphics = ~true;
    export_to_abaqus= true; SurfElType ='SFM3D4';;
    ElType ='C3D8';
    %             ElType ='C3D8H';
    ElType ='C3D8I';
    %     %     %     ElType ='C3D8IH';
%     ElType ='C3D8RH';
                                %                     ElType ='C3D20RH'; SurfElType ='SFM3D8';;
    
    Rm= rotmat(aangle/180*pi* [0,0,1]);
    
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
          
    eltyd(eix).description ='H8MSGSO';
    eltyd(eix).mf =@(ref)H8_block(a,b,t,ref*na,ref*nb,ref*nt);
    eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
        'Rm',Rm,'integration_rule',gauss_rule(struct('dim',3,'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eix=eix+1;
    
    %         eltyd(eix).description ='H27-BbarX';
    %         eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear_bbarx(struct('fes',fes, 'material',mater,...
    %             'Rm',Rm,...
    %             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',3)),...
    %             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
    %             'pv_bfun',@(p)[[1;p(1);p(2);p(3)],[0;0;0;0];[0;0;0;0],[1;p(1);p(2);p(3)]],'nconstrained',2));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
    %         eix=eix+1;
    %
    %     eltyd(eix).description ='H27-BbarX';
    %     eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbarx(struct('fes',fes, 'material',mater,...
    %         'Rm',Rm,...
    %         'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',3)),...
    %         'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
    %         'pv_bfun',@(p)[polynomial_basis(3,1,p),0*polynomial_basis(3,1,p);...
    %         0*polynomial_basis(3,1,p),polynomial_basis(3,1,p)],'nconstrained',2));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
    %     eix=eix+1;
    
    legends ={};
    for eix = 1:length(eltyd)
        ns=[];
        uzs=[];
        uys=[];
        energies=[];
        
        %         for ref = [4,8,16]
        for ref = refs
            
            % Mesh
            [fens,fes] = eltyd(eix).mf(ref);
                        count(fens)
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
            
            
            nc=[fenode_select(fens, struct('box', [a,a,0,0,0,0],...
                'inflate',tolerance))];
            ab_displacement=gather_values(model_data.u,nc);
            abz =mean(ab_displacement(:,3) );
            aby =mean(ab_displacement(:,2) );
            energy =model_data.work;
            
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
                AE.NSET_NSET('CORN',nc);
                AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
                AE.ORIENTATION('Rm', Rm(:,1), Rm(:,2));
                AE.SOLID_SECTION('Material','Rm','All','Hourglass');
                AE.HOURGLASS('STIFFNESS',5e6);
                %                 AE.SOLID_SECTION('Material','Rm','All');
                AE.SURFACE_SECTION('Traction');
                AE.NSET_NSET('xfix',find(model_data.u.is_fixed(:,1)));
                AE.NSET_NSET('yfix',find(model_data.u.is_fixed(:,2)));
                AE.NSET_NSET('zfix',find(model_data.u.is_fixed(:,3)));
                AE.END_INSTANCE();
                AE.END_ASSEMBLY();
                AE.MATERIAL('Material');
                AE.ELASTIC(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
                %               AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
                  AE.SECTION_CONTROLS('Hourglass','Hourglass = stiffness');
                AE.STEP_PERTURBATION('Linear solution');
                AE.DLOAD('ASSEM1.INSTNC1.Traction',model_data.boundary_conditions.traction{1}.traction);
                AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
                AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
                AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
                AE.NODE_PRINT('ASSEM1.INSTNC1.CORN');
                AE.ENERGY_PRINT();
                AE.END_STEP();
                AE.close();
                %                 delete([AE.filename '.dat']);
                system(['abaqus job=' [AE.filename ]]);
                pause(5);
                try
                    d =[0,0,0];
                    %                     d= extract_displacement_from_abaqus_dat([AE.filename '.dat'],'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_CORN');
                    energy= extract_energy_from_abaqus_dat([AE.filename '.dat'],'APPROXIMATE ENERGY TOTALS -');
                catch,
                    disp('failed');
                    d=[0,0, 0]; energy = 0;
                end
                abz=d(3);
                aby=d(2);
                eltyd(eix).description=ElType;
            end
            
            ns=[ns,model_data.u.nfreedofs];
            uzs=[uzs,abz];
            uys=[uys,aby];
            energies  = [energies,energy];
            
            if (graphics )
                figure
                model_data.postprocessing.u_scale= u_scale;
                model_data.postprocessing.outputRm=eye(3);
                model_data.postprocessing.cmap= abaquscolors;
                model_data.postprocessing.boundary_only= true;
                model_data.postprocessing.stress_component=6;
                model_data.postprocessing.stress_range=[-600, 600];
                model_data.postprocessing.camera= [ -0.551266963760385  -0.491703063433634   0.395765903053231   0.046327725040208   0.050000000000000   0.004472559585213   0.329913457292308   0.285915160851954 0.899671957711497   8.330593849827762];
                model_data=deformation_plot_stress_elementwise(model_data);
                axis tight
                %                 title(eltyd(eix).description)
                %                 model_data=deformation_plot_deformation(model_data);
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
         disp(num2str(energies,'%15.12e '));
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
    if (~graphics )
        legend (legends);
        title (['Fiber dir =[' num2str(Rm(:,1)') ']'])
        labels(  'Number of equations', 'Estimated true error')
        grid on
        
        set_decades_on_axis (gca)
        set_pub_defaults
    end
end
