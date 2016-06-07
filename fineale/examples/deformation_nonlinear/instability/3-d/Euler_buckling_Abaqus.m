%
function Euler_buckling
E=100;
nu=0.3;
L= 28; % Length of the beam
W = 1; % width
H = 3; % height
% Uniform compressive traction load for fundamental buckling mode
magn_cr = E*min([H,W])^3*max([H,W])/12*pi^2/(2*L)^2/W/H;
magn=magn_cr;
magn_cr_other = E*max([H,W])^3*min([H,W])/12*pi^2/(2*L)^2/W/H;
nW=4; nH=4;
graphics= ~true;  scale=5;
export_to_abaqus= ~true; SurfElType ='SFM3D4';;
ElType ='C3D8';
    ElType ='C3D8H';
% ElType ='C3D8I';
%     ElType ='C3D8IH';
    ElType ='C3D8RH';
%         ElType ='C3D20R'; SurfElType ='SFM3D8';;
% ElType ='C3D20R'; SurfElType ='SFM3D8';;

rand('state',0);% try to comment out this line and compare
%                   results for several subsequent runs
%
    

eix=1;
clear eltyd



eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).blf =@femm_deformation_nonlinear_h8msgso;
eltyd(eix).integration_rule=gauss_rule(struct('dim',3,'order',2));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
eltyd(eix).styl='ko-';
eix=eix+1;

% eltyd(eix).description='H8';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3,'order',2));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% eltyd(eix).styl='ko-';
% eix=eix+1;

% eltyd(eix).description ='H64';
% eltyd(eix).mf =@h8H64;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',4));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
% eltyd(eix).styl='r*-';
% eix=eix+1;
%
% eltyd(eix).description ='H27';
% eltyd(eix).mf =@h8H27;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',3));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% eltyd(eix).styl='ks-';
% eix=eix+1;
% eltyd(eix).description ='H20R';
% eltyd(eix).mf =@H20_block;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',2));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% eltyd(eix).styl='ks-';
% eix=eix+1;

    function [Buckle] =Compute
        %          Mesh
        [fens,fes] = eltyd(eix).mf (W, H, L, nW, nH, nL);
        % Material
        prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',1.0));
        mater = material_deformation_stvk_triax (struct('property',prop ));
        femm = eltyd(eix).blf (struct ('material',mater, 'fes',fes,'integration_rule',eltyd(eix).integration_rule));
        bfes =mesh_boundary(fes, []);
        icl = fe_select(fens, bfes, struct('box', [-Inf,Inf,-Inf,Inf,L,L],'inflate',1/1000));
        %          gv=drawmesh( {fens,fes},'facecolor','none'); gv=drawmesh( {fens,bfes(icl)},'gv',gv,'facecolor','blue');
        %         view(3); pause(1)
        % Finite element block
        femm = eltyd(eix).blf (struct ('material',mater, 'fes',fes,...
            'integration_rule',eltyd(eix).integration_rule));
        efemm = eltyd(eix).blf (struct ('material',mater, 'fes',subset(bfes,icl),...
            'integration_rule',eltyd(eix).surface_integration_rule));
        sfemm = eltyd(eix).blf (struct ('material',mater, 'fes',bfes,...
            'integration_rule',eltyd(eix).surface_integration_rule));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = 0*geom; % zero out
        % Apply EBC's
        %  Clamped end
        ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0,],'inflate',1/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=[];
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        
        u   = apply_ebc (u);
        
        % Number equations
        u   = numberdofs (u);
        u = u*0; % zero out the displacement
        
       
        femm  =associate_geometry(femm,geom); 
        femm1=femm;                                           % Make a copy of the state
        fi=force_intensity(struct('magn', [0;0; -magn]));
        F =  distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
        %         sum (get(F,'vec'))
        
        % Update the FEMM
        K = stiffness(femm1, sysmat_assembler_sparse, geom, u,0*u);
        u = scatter_sysvec(u, K\F); % Displacement increment
        lnl=fenode_select(fens,struct('box',[-inf inf 0 0 -inf inf]));
        load_displacement=mean(gather_values(u,lnl));
        if (graphics)
            gv=graphic_viewer;
            gv=reset (gv,[]);
            draw(femm,gv, struct ('x', geom,'u',scale*u));
            draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
        end
        
        Kgeo = stiffness_geo(femm1, sysmat_assembler_sparse, geom, u,0*u);
        
        neigvs = 4;
        options.tol =1e-206;
        options.maxit = 5000;
        options.disp = 0;
        [Phi,Omega]=eigs(-(Kgeo+Kgeo')/2,(K+K')/2,...
            neigvs,'lm', options);
        %         diag(Omega)
        %         [Omegas,ix]=sort(1./diag(Omega)) ;
        Omegas=(1./(diag(Omega)));
        Buckle =Omegas;
        
        if (export_to_abaqus)
            date_now = clock;
            s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
            AE=Abaqus_exporter;
            AE.open([mfilename '-' ElType '-' num2str(nL) '-' s '.inp']);
            AE.HEADING([mfilename ' ' 'ElType=' ElType ' ' 'nL=' num2str(nL)]);
            AE.PART('PART1');
            AE.END_PART();
            AE.ASSEMBLY('ASSEM1');
            AE.INSTANCE('INSTNC1','PART1');
            AE.NODE(geom.values);
            AE.ELEMENT(ElType,'All',1,femm1.fes.conn);
            AE.ELEMENT(SurfElType,'Traction',count(femm1.fes)+1,efemm.fes.conn);
            AE.ORIENTATION('Global', [1,0,0], [0,1,0]);
            AE.SOLID_SECTION('Material','Global','All','Hourglass');
            %             AE.SOLID_SECTION('Material','Global','All');
            AE.SURFACE_SECTION('Traction');
            AE.NSET_NSET('xfix',find(u.is_fixed(:,1)));
            AE.NSET_NSET('yfix',find(u.is_fixed(:,2)));
            AE.NSET_NSET('zfix',find(u.is_fixed(:,3)));
            AE.END_INSTANCE();
            AE.END_ASSEMBLY();
            AE.MATERIAL('Material');
            AE.ELASTIC_ISOTROPIC(E,nu);
            %             AE.SECTION_CONTROLS('Hourglass','Hourglass = enhanced');
            AE.STEP_PERTURBATION_BUCKLE('buckle',neigvs);
            AE.DLOAD('ASSEM1.INSTNC1.Traction',fi.magn);
            AE.BOUNDARY('ASSEM1.INSTNC1.xfix',1);
            AE.BOUNDARY('ASSEM1.INSTNC1.yfix',2);
            AE.BOUNDARY('ASSEM1.INSTNC1.zfix',3);
            AE.END_STEP();
            AE.close();
            %                 delete([AE.filename '.dat']);
            abaqus_job([AE.filename ]);
            AW=Abaqus_lck_watcher();
            AW.wait([AE.filename '.lck']);
            try
                AR=Abaqus_dat_reader();
                d= AR.extract_buckling_from_abaqus_dat([AE.filename '.dat'],...
                    'MODE NO      EIGENVALUE',neigvs)
            catch,
                d=zeros(1,neigvs); 
            end
            Buckle =d;
        end
        
        for i=1:neigvs
            disp(['  Eigenvector ' num2str(i) ' eigenvalue ' num2str(Buckle(i)) ]);
        end
        
        % Plot
        if graphics
            scale=100;
            for i=(1:neigvs)
                phi = scatter_sysvec(u,Phi(:,i));
                phiv = phi.values;
                %             norm(K*gather_sysvec(phi) +Omegas(i)*Kgeo*gather_sysvec(phi))
                %             disp(['  Buckling load '  num2str((Omegas(i)*Fmag)) ]);
                %     clf;
                if (graphics)
                    gv=graphic_viewer;
                    gv=reset (gv,[]);
                    scale =L/3/max(abs(gather_sysvec(phi)));
                    phivmag = sqrt(phiv(:,1).^2+phiv(:,2).^2+phiv(:,3).^2);
                    dcm=data_colormap(struct ('range',[min(phivmag),max(phivmag)], 'colormap',colormap('parula')));
                    colors=map_data(dcm, phivmag);
                    colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
                    %                 set(gca,'FontSize',16)
                    %     camset(gv,[-0.1591    0.0085   -0.3252    0.0001    0.0043    0.0013   -0.8985    0.0237    0.4383    5.4322]);
                    draw(femm,gv, struct ('x', geom,'u',0*phi, 'facecolor','none'));
                    draw(femm,gv, struct ('x', geom,'u',+scale*phi,'colorfield',colorfield));
                    %     text(0,1.1*R,1.1*R,['\omega_' num2str(i) '=' num2str(sqrt(Omega(i,i))/2/pi) ],'FontSize',24);
                    grid on; set_graphics_defaults
                                    pause(1)
                end
            end
        end
    end

for eix = 1:length(eltyd)
    nW=1; nH=2;
    nelperedge =(2:2:8)+12;
    nelperedge = 20;
    Buckles=[];
    for nL =  nelperedge
        [Buckle] = Compute;
        Buckles=[Buckles,Buckle];
    end
    disp(['Data{end+1}=[']);
    disp([num2str(nelperedge,15)]);
    disp([num2str(Buckles,15)]);
    if (export_to_abaqus)
        disp(['];' 'Style{end+1}=''' ElType ''';' 'Legends{end+1} =''' ElType ''';']);
    else
        disp(['];' 'Style{end+1}=''' name_to_style(eltyd(eix).description) ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
    end
    
    
end
end