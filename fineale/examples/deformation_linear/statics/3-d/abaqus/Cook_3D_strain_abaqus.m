function  Cook_3D_strain
% Parameters:
lambda = 0.75e4;
G= 0.375;
convutip=16.432584376299037;% NICE-T6 t3block2d with Richardson
nu =1/2*lambda/(G+lambda);
E=G*(3*lambda+2*G)/(G+lambda);
tol=48/100000;
thickness=48/48*12;% Thick sample
thickness=48/48*0.1;% Thin sample
magn=1;
aspect=1;
graphics = ~true;
u_scale=1;
export_to_abaqus= true; SurfElType ='SFM3D4';;
ElType ='C3D8';
ElType ='C3D8H';
ElType ='C3D8I';
ElType ='C3D8IH';
%     ElType ='C3D8';
%          ElType ='C3D8RH';
%     ElType ='C3D20R'; SurfElType ='SFM3D8';;

ns= [2,4,8,16,32,64,128];
ns=[2,4,8,16];
ns= [2,4,6];


format long
% clc
disp(['% ' mfilename]);
disp(['nu=' num2str(nu,12) '; ']);

prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;

eltyd(eix).description ='H8';% tetrahedron
eltyd(eix).mf =@H8_block;
eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

% eltyd(eix).description ='H8MSGSO';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
%     'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

% eltyd(eix).description ='H8MSGS';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,'material',mater,...
%     'integration_rule', gauss_rule(struct('dim',3','order',2)), 'stab_fraction', 0.1,'weights',(1-2*nu)/2));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

% eltyd(eix).description ='H8-3FS';% tetrahedron
% eltyd(eix).mf =@H8_block;
% O=@(p)polynomial_basis(3,0,p);
% eltyd(eix).femmf =@(fes)femm_deformation_linear_3fs(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'pv_bfun',@(p)[O(p),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1);
%     zeros(1,1),O(p),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1);
%     zeros(1,1),zeros(1,1),O(p),zeros(1,1),zeros(1,1),zeros(1,1);
%     zeros(1,1),zeros(1,1),zeros(1,1),O(p),zeros(1,1),zeros(1,1);
%     zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),O(p),zeros(1,1);
%     zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),O(p);],...
%     'nconstrained',6,'psi',[0.2,0.2,0.2, [1,1,1]]/2.5));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

% eltyd(eix).description ='H20R';
% eltyd(eix).mf =@H20_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%     'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
% %
%
% eltyd(eix).description ='H20-Bbar';
% eltyd(eix).mf =@H20_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
%    'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
%     'pv_bfun',@(p)[1;p(1);p(2);p(3)],'nconstrained',1));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
%     eltyd(eix).description ='T10';% tetrahedron
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
%         'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;
%
% % eltyd(eix).description ='THEX';% tetrahedron
% % eltyd(eix).mf =@T15_blocka;
% % eltyd(eix).femmf =@(fes)femm_deformation_linear_thex(struct('fes',fes,'material',mater,...
% %     'stabfact',0.01*E));
% % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% % eix=eix+1;
%
% eltyd(eix).description ='H8-Bbar';% tetrahedron
% eltyd(eix).mf =@H8_block;plate_2dir_strong_abaqus_comparison
% eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),'nconstrained',1));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H8-BbarX';% tetrahedron
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_bbarx(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'pv_bfun',@(p)[polynomial_basis(3,0,p),0,0;    0,polynomial_basis(3,0,p),0; 0,0,polynomial_basis(3,0,p)],'nconstrained',3,'psi',1-0.5+nu));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% % eltyd(eix).description ='H27-Bbar';% tetrahedron
% % eltyd(eix).mf =@H27_block;
% % eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
% %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
% %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
% %     'pv_bfun',@(p)[1;p(1);p(2);p(3);]));%@(p)[1;p(1);p(2);p(3);]
% % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% % eix=eix+1;
%
%     eltyd(eix).description ='T10-Bbar';
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
%         'integration_rule_constrained',tet_rule(struct('npts',1)),...
%         'integration_rule_unconstrained',tet_rule(struct('npts',4)),...
%          'pv_bfun',@(p)[1]));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;
%
% %         eltyd(eix).description ='T10-Bbar';
% %     eltyd(eix).mf =@T10_block;
% %     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
% %         'integration_rule_constrained',tet_rule(struct('npts',1)),...
% %         'integration_rule_unconstrained',tet_rule(struct('npts',4)),...
% %          'pv_bfun',@(p)[1]));
% %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% %     eix=eix+1;
%
%
% % eltyd(eix).description ='NICE-H8';% tetrahedron
% % eltyd(eix).mf =@H8_block;
% % eltyd(eix).femmf =@(fes)femm_deformation_linear_nice(struct('fes',fes,'material',mater,...
% %     'integration_rule',tensprod_nq_rule(struct('dim',3, 'order',1))));
% % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% % eix=eix+1;
%
% % eltyd(eix).description ='H8-GSRI';% tetrahedron
% % eltyd(eix).mf =@H8_block;
% % eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
% %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
% %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% % eix=eix+1;
% %
% % eltyd(eix).description ='H8-SRI';% tetrahedron
% % eltyd(eix).mf =@H8_block;
% % eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes,'material',mater,...
% %     'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
% %     'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
% % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% % eix=eix+1;
%
%
% %     eltyd(eix).description ='T10';% tetrahedron
% %     eltyd(eix).mf =@T10_block;
% %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
% %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% %     eltyd(eix).styl='k^--';
% %     eix=eix+1;
% % %
% %     eltyd(eix).description ='T10-GSRI';
% %     eltyd(eix).mf =@T10_block;
% %     eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
% %         'integration_rule_constrained',tet_rule(struct('npts',1)),...
% %         'integration_rule_unconstrained',tet_rule(struct('npts',4))));
% %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% %     eix=eix+1;
%
% %     eltyd(eix).description ='NICE-T10-GSRI';% tetrahedron
% %     eltyd(eix).mf =@T10_block;
% %     eltyd(eix).femmf =@(fes)fem_deformation_linear_nice_gsri(struct('fes',fes,'material',mater,...
% %         'integration_rule_constrained',simplex_nq_rule(struct('dim',3, 'order',2)),...
% %         'integration_rule_unconstrained',tet_rule(struct('npts',4))));
% %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% %     eltyd(eix).styl='c^-';
% %     eix=eix+1;
%
%
% %         eltyd(eix).description ='H27-GSRI';% tetrahedron
% %         eltyd(eix).mf =@H27_block;
% %         eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
% %             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
% %             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3))));
% %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% %         eix=eix+1;
% %
% %
% %         eltyd(eix).description ='H64-GSRI';
% %         eltyd(eix).mf =@H64_block;
% %         eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
% %             'integration_rule_constrained', hex_27pt_rule(),...
% %             'integration_rule_unconstrained',hex_64pt_rule()...
% %             ));
% %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
% %         eltyd(eix).styl='kx-';
% %         eix=eix+1;
%
% %     eltyd(eix).description ='H20-GSRI-8-9g';
% %     eltyd(eix).mf =@H20_block;
% %     eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
% %         'integration_rule_constrained', hex_8pt_rule(),...
% %         'integration_rule_unconstrained',hex_9ptg_rule()...
% %         ));
% %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
% %     eltyd(eix).styl='bd-';
% %     eix=eix+1;
% %
%
% %
% %     eltyd(eix).description ='H20';
% %     eltyd(eix).mf =@H20_block;
% %     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
% %         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
% %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% %     eltyd(eix).styl='r.--';
% %     eix=eix+1;
% %
% %     eltyd(eix).description ='H27';
% %     eltyd(eix).mf =@H27_block;
% %     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
% %         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
% %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% %     eltyd(eix).styl='rh--';
% %     eix=eix+1;
% %
% %         eltyd(eix).description ='H20-GSRI-8-21';
% %         eltyd(eix).mf =@H20_block;
% %         eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
% %             'integration_rule_constrained', hex_8pt_rule(),...
% %             'integration_rule_unconstrained',hex_21pt_rule()...
% %             ));
% %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
% %         eltyd(eix).styl='kx-';
% %         eix=eix+1;
% %
% %         eltyd(eix).description ='H27-GSRI-8-27';
% %         eltyd(eix).mf =@H27_block;
% %         eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
% %             'integration_rule_constrained', hex_8pt_rule(),...
% %             'integration_rule_unconstrained',hex_27pt_rule()...
% %             ));
% %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
% %         eltyd(eix).styl='kd-';
% %         eix=eix+1;
% %
% %     eltyd(eix).description ='H20R-9';
% %     eltyd(eix).mf =@H20_block;
% %     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
% %         'integration_rule',hex_9pt_rule()));
% %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% %     eltyd(eix).styl='ks-';
% %     eix=eix+1;
% %


legends ={};
for eix = 1:length(eltyd)
    uzs=[];
    ndofs = [];
    
    for n=ns  %,16,24,32,48
        % Create the mesh and initialize the geometry
        [fens,fes]= eltyd(eix).mf(48,44,thickness, n, aspect*n,1);%max([round(nt*ref/2),2]));
        dxy=min(48,44)/n/aspect/100;
        xy=fens.xyz;
        for i=1:count (fens)
            xy(i,:)=[xy(i,1),xy(i,2) + (xy(i,1)/48)*(44 -xy(i,2)/44*28),xy(i,3)];
        end
        fens.xyz=xy;
        
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
        essential.node_list = fenode_select (fens,struct ('box',[0,0,0, 44,0,thickness],'inflate',tol));
        model_data.boundary_conditions.essential{1} = essential;
        
        clear essential
        essential.component= [3];
        essential.fixed_value= 0;
        essential.node_list = [fenode_select(fens,struct('box',[-inf,inf,-inf,inf,0,thickness],'inflate',tol))];
        model_data.boundary_conditions.essential{2} = essential;
        
        clear traction
        bdry_fes = mesh_boundary(fes, []);
        bcl = fe_select(fens, bdry_fes, ...
            struct ('box',[48,48,44,60,0,thickness],'inflate',tol));
        traction.fes =subset(bdry_fes,bcl);;
        traction.traction= [0;magn/(60-44);0];
        traction.integration_rule =eltyd(eix).surface_integration_rule;
        model_data.boundary_conditions.traction{1} = traction;
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        if (graphics )
            model_data.postprocessing.u_scale= u_scale;
            %             model_data=deformation_plot_deformation(model_data);
            model_data.postprocessing.component= 'pressure';
            model_data.postprocessing.stress_range = [-0.3, 0.13];
            %              model_data=deformation_plot_stress(model_data);
            model_data=deformation_plot_stress_elementwise(model_data);
            %              title (eltyd(eix).description)
            axis tight
            view(2);
            axis tight
        end
        
        enl=fenode_select (fens,struct ('box',[48,48,52,52,0,thickness],'inflate',tol));
        uv=gather_values (model_data.u,enl);
        utip=mean(uv(:,2));
        
        
        if (export_to_abaqus)
            date_now = clock;
            s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
            AE=Abaqus_exporter;
            AE.open([mfilename '-' ElType '-' num2str(n) '-' s '.inp']);
            AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'n=' num2str(n)]);
            AE.PART('part1');
            AE.END_PART();
            AE.ASSEMBLY('ASSEM1');
            AE.INSTANCE('INSTNC1','PART1');
            AE.NODE(model_data.geom.values);
            AE.ELEMENT(ElType,'All',1,fes.conn);
            AE.ELEMENT(SurfElType,'Traction',count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
            AE.NSET_NSET('Corner',enl);
            AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
            AE.SOLID_SECTION('Material','Global','All','Hourglass');
            AE.SURFACE_SECTION('Traction');
            AE.NSET_NSET('xfix',find(model_data.u.is_fixed(:,1)));
            AE.NSET_NSET('yfix',find(model_data.u.is_fixed(:,2)));
            AE.NSET_NSET('zfix',find(model_data.u.is_fixed(:,3)));
            AE.END_INSTANCE();
            AE.END_ASSEMBLY();
            AE.MATERIAL('Material');
            AE.ELASTIC_ISOTROPIC(E,nu);
            AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
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
                    'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_CORNER',...
                    length(enl));
            catch,
                d=[0,0, 0]; energy = 0;
            end
            utip=mean(d(:,2));
            eltyd(eix).description=ElType;
        end
        
        disp (['%     Normalized displacement =' num2str(utip/convutip)])
        uzs =[uzs utip] ;
        ndofs = [ndofs,model_data.u.nfreedofs];
    end
    legends{end+1} =eltyd(eix).description;
    format long
    disp(['Data{end+1}=[']);
    disp([num2str(ns)]);
    disp([num2str(ndofs)]);
    disp([num2str(uzs,12)]);
    disp(['];' ''''  ''';' 'Description{end+1} =''' eltyd(eix).description ''';']);
    if (~graphics )
        semilogx(ndofs,uzs/convutip,name_to_style(eltyd(eix).description),'linewidth',3);  hold on
    end
end
if (~graphics )
    legend (legends);
    grid on
    labels(  'Number of degrees of freedom', 'Normalized deflection')
else
    camset(model_data.postprocessing.gv,1.0e+02 *[-0.845458458438780  -2.592911972715104   3.139339064298713   0.223225118317536   0.385368318712585  0.005000000000000   0.001127144822962   0.007008300749704   0.007043668444054   0.079654194001118]);
end
set_graphics_defaults
end