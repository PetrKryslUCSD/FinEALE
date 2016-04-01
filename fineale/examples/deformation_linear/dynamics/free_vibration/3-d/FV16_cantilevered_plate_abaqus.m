% This is a test recommended by the National Agency for Finite Element Methods and Standards (U.K.):
% Test FV16 from NAFEMS publication TNSB, Rev. 3, “The Standard NAFEMS Benchmarks,” October
% 1990.
% Model: Plate thickness = 0.05 m.
% Material: Young’s modulus = 200 GPa, Poisson’s ratio = 0.3, density = 8000 kg/m 3 .
% Boundary conditions: along the y-axis.
% NAFEMS 0.421 1.029 2.582 3.306 3.753 6.555
function FV16_cantilevered_plate_abaqus
pu=physical_units_struct;
% Parameters:
E = 200e9*pu.PA;
nu = 0.3;
rho= 8000*pu.KG/pu.M^3;
Rm=eye(3);
a=10*pu.M; b=a; h= 0.05*pu.M;
%     Number of eigenvalues
neigvs=6;
% Mesh refinement
n1=10;
na= n1; nb= n1; nh =4;
htol=h/nh/1000;
export_to_abaqus=   true; SurfElType ='SFM3D4';;
    ElType ='C3D8';
    %             ElType ='C3D8H';
                        ElType ='C3D8I';
%    ElType ='C3D8IH';
% ElType ='C3D8R';
    %                 ElType ='C3D20R'; SurfElType ='SFM3D8';;
    
% Fundamental frequency
graphics = ~false;

prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
mater = material_deformation_linear_triax (struct('property',prop ));

% prop_stab=property_deformation_linear_iso(struct('E',stab_fraction*E/2/(1+nu),'nu',0,'rho',rho));
% mater_stab = material_deformation_linear_triax (struct('property',prop_stab ));

% B0=@(p)polynomial_basis(3,0,p);
% switch nconstrained
%     case 3
%         PV=@(p)[B0(p),0*B0(p),0*B0(p); 0*B0(p),B0(p),0*B0(p); 0*B0(p),0*B0(p),B0(p)];
%     case 2
%         PV=@(p)[B0(p),0*B0(p); 0*B0(p),B0(p)];
%     case 1
%         PV=@(p)[B0(p)];
%     otherwise
%         error
% end

clear eltyd
eix=1;

   
eltyd(eix).description ='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

% eltyd(eix).description ='H20-Bbar';
%     eltyd(eix).mf =@H8_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
%         'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%         'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2)),...
%         'pv_bfun',PV,'nconstrained',nconstrained));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eix=eix+1;

% B1=@(p)polynomial_basis(3,1,p); PV=@(p)[B1(p)];
%     eltyd(eix).description ='H27-Bbar';
%     eltyd(eix).mf =@H27_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
%         'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',3)),...
%         'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',3)),...
%         'pv_bfun',PV,'nconstrained',nconstrained));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eix=eix+1;
%
% eltyd(eix).description ='H27-Bbar';% tetrahedron
% eltyd(eix).mf =@H27_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',3)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
%     'integration_rule',gauss_rule(struct('dim',3, 'order',3)),...
%         'pv_bfun',@(p)[1;p(1);p(2);p(3);]));%@(p)[1;p(1);p(2);p(3);]
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% eix=eix+1;

% %     % Selective reduced integration hexahedron
%     eltyd(eix).description ='H8-SRI';
%     eltyd(eix).mf =@H8_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
%         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2)),...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%     eix=eix+1;
%
%     eltyd(eix).description ='H20R';
%     eltyd(eix).mf =@H20_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eix=eix+1;


legends ={};
for eix = 1:length(eltyd)
    
    ns=[];
    uzs=[];
    
    %% Create the mesh and initialize the geometry
    [fens,fes]= eltyd(eix).mf(a,b,h,na,nb,nh);
    %                     drawmesh({fens,fes},'shrink', 0.8,'facecolor','red');
    
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
    essential.node_list = fenode_select(fens,struct('box',[0 0 0 b 0 h],'inflate',htol));
    model_data.boundary_conditions.essential{1} = essential;
    
    model_data.neigvs= neigvs;
    model_data.omega_shift=(2*pi*0.01) ^ 2;
    model_data.use_factorization= true;
    model_data.use_lumped_mass=~true;
    % Solve
    model_data = deformation_linear_modal_analysis(model_data);
    
    if (export_to_abaqus)
        date_now = clock;
        s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
        AE=Abaqus_exporter;
        AE.open([mfilename '-' ElType '-' num2str(n1) '-' s '.inp']);
        AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'n1=' num2str(n1)]);
        AE.PART('part1');
        AE.END_PART();
        AE.ASSEMBLY('ASSEM1');
        AE.INSTANCE('INSTNC1','PART1');
        AE.NODE(model_data.geom.values);
        AE.ELEMENT(ElType,'All',1,fes.conn);
        %         AE.ELEMENT(SurfElType,'Traction',count(fes)+1,model_data.boundary_conditions.traction{1}.fes.conn);
        %         AE.NSET_NSET('MIDEDGE',nc);
        AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
        AE.ORIENTATION('Rm', Rm(:,1), Rm(:,2));
        AE.SOLID_SECTION('Material','Rm','All','Hourglass');
        %                 AE.HOURGLASS('STIFFNESS',0.3/10);
        %                 AE.SOLID_SECTION('Material','Rm','All');
        AE.SURFACE_SECTION('Traction');
        AE.NSET_NSET('xfix',find(model_data.u.is_fixed(:,1)));
        AE.NSET_NSET('yfix',find(model_data.u.is_fixed(:,2)));
        AE.NSET_NSET('zfix',find(model_data.u.is_fixed(:,3)));
        AE.END_INSTANCE();
        AE.END_ASSEMBLY();
        AE.MATERIAL('Material');
        AE.ELASTIC_ISOTROPIC(E,nu);
        AE.DENSITY(rho);
        %                     AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
          %    AE.SECTION_CONTROLS('Hourglass','Hourglass = stiffness');
        AE.STEP_FREQUENCY(neigvs);
        %         AE.DLOAD('ASSEM1.INSTNC1.Traction',model_data.boundary_conditions.traction{1}.traction);
        AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
        AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
        AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
        %         AE.NODE_PRINT('ASSEM1.INSTNC1.MIDEDGE');
        AE.ENERGY_PRINT();
        AE.END_STEP();
        AE.close();
        %                 delete([AE.filename '.dat']);
        abaqus_job([AE.filename ]);
        AW=Abaqus_lck_watcher();
        AW.wait([AE.filename '.lck']);
        try
            d= extract_displacement_from_abaqus_dat([AE.filename '.dat'],'THE FOLLOWING TABLE IS PRINTED FOR NODES BELONGING TO NODE SET ASSEM1_INSTNC1_MIDEDGE');
        catch,
            d=[0,0, 0]; energy = 0;
        end
        abz=d(3);
        aby=d(2);
        eltyd(eix).description=ElType;
    end
    
    f=model_data.Omega'/2/pi;
    disp ([eltyd(eix).description  ': '])
    disp(['Frequency ' num2str(abs(f)) ' [Hz]' ]);
    disp(['Frequency ' num2str((f)) ' [Hz]' ]);
    disp(['             f/f_analytical=' num2str(f./[0.421 1.029 2.582 3.306 3.753 6.555]*100) '%']);
    %         clf;
    if (graphics)
        model_data.postprocessing.u_scale= 2;
        model_data.postprocessing.colormap= gray;cadcolors2;
        model_data.postprocessing.modelist= 1:model_data.neigvs;
        model_data.postprocessing.save_frame= true;
        model_data.postprocessing.frame_name= [eltyd(eix).description '-n1-' num2str(n1)];
        model_data=deformation_plot_modes(model_data);
    else
        plot(model_data.Omega/(2*pi),name_to_style(eltyd(eix).description));;  hold on
        figure(gcf);
    end
    legends{end+1} =eltyd(eix).description;
end
if (graphics)
else
    legend (legends);
    labels(  'Mode number', 'Natural frequency')
    grid on
    set_graphics_defaults
end
end
