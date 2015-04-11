function  Twisted_beam_abaqus
% Parameters:
E=0.29e8;
nu=0.22;
W=1.1;
L=12;
t= 0.32;
nl=2; nt=1; nw=1;
p=  1/W/t;
loadv=[0;0;p];dir=3;uzex=0.005424534868469; % Harder: 5.424e-3;
loadv=[0;p;0];dir=2;uzex=0.001753248285256; % Harder: 1.754e-3;
graphics = true;
u_scale=1000;
export_to_abaqus= true; SurfElType ='SFM3D4';;
    ElType ='C3D8';
    ElType ='C3D8H';
                ElType ='C3D8I';
    %     % ElType ='C3D8IH';
                ElType ='C3D8RH';
    %     ElType ='C3D20R'; SurfElType ='SFM3D8';;

prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

prop_stab=property_deformation_linear_iso(struct('E',0.1*E/2/(1+nu),'nu',0));
mater_stab = material_deformation_linear_triax (struct('property',prop_stab ));

clear eltyd
eix=1;

%               eltyd(eix).description ='H8MSGS';
%               eltyd(eix).mf =@H8_block;
%         eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,'material',mater,'nconstrained',1, 'stab_fraction',0.1));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%         eix=eix+1;

eltyd(eix).description ='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

%           eltyd(eix).description ='H8MSGS';
%     eltyd(eix).mf =@H8_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2)),...
%         'material',mater, 'stab_fraction',0.1,'weights',[]));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eix=eix+1;
    

    %           eltyd(eix).description ='H8MSGS';
    %     eltyd(eix).mf =@H8_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgs(struct('fes',fes,...
    %         'stabilization_material',mater_stab,'material',mater));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eix=eix+1;
    

% eltyd(eix).description ='THEX';
% eltyd(eix).mf =@T15_blocka;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_thex(struct('fes',fes,'material',mater,...
%     'stabfact',0.01*E));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%


%         eltyd(eix).description ='H8-3FS';
%         eltyd(eix).mf =@H8_block;
%         O=@(p)polynomial_basis(3,0,p);
%         eltyd(eix).femmf =@(fes)femm_deformation_linear_3fs(struct('fes',fes,'material',mater,...
%             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
%             'pv_bfun',@(p)[O(p),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1);
%             zeros(1,1),O(p),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1);
%             zeros(1,1),zeros(1,1),O(p),zeros(1,1),zeros(1,1),zeros(1,1);
%             zeros(1,1),zeros(1,1),zeros(1,1),O(p),zeros(1,1),zeros(1,1);
%             zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),O(p),zeros(1,1);
%             zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),zeros(1,1),O(p);],...
%             'nconstrained',6,'psi',[0.2,0.2,0.2, [1,1,1]]/20));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%         eix=eix+1;
    
% eltyd(eix).description ='NICE-H8';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_nice(struct('fes',fes,'material',mater,...
%     'integration_rule',tensprod_nq_rule(struct('dim',3, 'order',1))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H8-GSRI';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H20R';
% eltyd(eix).mf =@H20_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%     'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

%             eltyd(eix).description ='TMH27';
%             eltyd(eix).mf =@H27_block;
%             eltyd(eix).femmf =@(fes)fem_deformation_linear_tmh27(struct('fes',fes, 'material',mater,...
%             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3))));
%             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%             eltyd(eix).styl='kp-';
%             eix=eix+1;
%
%             eltyd(eix).description ='H27-GSRI';
%             eltyd(eix).mf =@H27_block;
%             eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%             'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%             'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3))));
%             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%             eltyd(eix).styl='kd--';
%             eix=eix+1;
%
%     eltyd(eix).description ='H64';
%     eltyd(eix).mf =@H64_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%     eltyd(eix).styl='k*--';
%     eix=eix+1;
%
%     eltyd(eix).description ='H27';
%     eltyd(eix).mf =@H27_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%     eltyd(eix).styl='ks--';
%     eix=eix+1;
%
%     eltyd(eix).description ='CH27-NICE-GSRI';
%     eltyd(eix).mf =@H27_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear_ch27_nice_gsri(struct('fes',fes, 'material',mater));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='rs-';
%     eix=eix+1;
%
%     eltyd(eix).description ='T10';
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;
%
%     eltyd(eix).description ='T10-GSRI';
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%         'integration_rule_constrained',tet_rule(struct('npts',1)),...
%         'integration_rule_unconstrained',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eltyd(eix).styl='kv-';
%     eix=eix+1;
%
%     eltyd(eix).description ='H20R';
%     eltyd(eix).mf =@H20_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='ro--';
%     eix=eix+1;
%
%     %         % Selective reduced integration hexahedron
%     eltyd(eix).description ='H8-SRI';
%     eltyd(eix).mf =@H8_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
%         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='md--';
%     eix=eix+1;
%
mix = 1;
mesd(mix).ref=1:4;[1,2,4];
mix = mix+1;


for eix = 1:length(eltyd)
    ns=[];
    
    for mix = 1:length(mesd)
        uzs=[];
        neqns = [];
        
        for     ref=mesd(mix).ref;
            %% Create the mesh and initialize the geometry
            [fens,fes]= eltyd(eix).mf(L,W,t, nl*ref,nw*ref,nt*ref);%max([round(nt*ref/2),2]));
            xy=fens.xyz;
            for i=1:count (fens)
                a=xy(i,1)/L*(pi/2); y=xy(i,2)-(W/2); z=xy(i,3)-(t/2);
                xy(i,:)=[xy(i,1),y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
                %                         get(fens(i),'id'),get(fens(i),'xyz')
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
            essential.node_list = fenode_select (fens,struct ('box',[0 0 -100*W 100*W -100*W 100*W],'inflate',0.001*t));
            model_data.boundary_conditions.essential{1} = essential;
            
            clear traction
            bdry_fes = mesh_boundary(fes, []);
            bcl = fe_select(fens, bdry_fes, ...
                struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.0001*t));
            traction.fes =subset(bdry_fes,bcl);;
            traction.traction= loadv;
            traction.integration_rule =eltyd(eix).surface_integration_rule;
            model_data.boundary_conditions.traction{1} = traction;
            
            
            enl=fenode_select (fens,struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.01*t));
            
            % Solve
            model_data =deformation_linear_statics(model_data);
            uv=gather_values (model_data.u,enl);
            uz=mean(uv(:,dir));
            
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
% AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
AE.SECTION_CONTROLS('Hourglass','Hourglass = stiffness');
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
                pause(5);
                try
                    d= extract_displacement_from_abaqus_dat([AE.filename '.dat'],...
                        'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_CORNER',...
                        length(enl));
                catch,
                    d=[0,0, 0]; energy = 0;
                end
                uz=mean(d(:,dir));
                eltyd(eix).description=ElType;
            end
            
            if (graphics )
                model_data.postprocessing.u_scale= u_scale;
                model_data.postprocessing.show_mesh= 1;
                model_data=deformation_plot_deformation(model_data);
            end
            
            disp (['%     Normalized displacement =' num2str(uz/uzex)])
            uzs =[uzs uz] ;
            neqns = [neqns,model_data.u.nfreedofs];
            ns=[ns,nl*ref] ;
        end
    end
    format long
    disp(['Data{end+1}=[']);
    disp([ns;uzs]);
    disp(['];' ''''  ''';' 'Description{end+1} =''' eltyd(eix).description ''';']);
    plot(ns,uzs,name_to_style(eltyd(eix).description));  hold on
end
end
