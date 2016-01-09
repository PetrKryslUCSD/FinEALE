function distorted
% Numerical inf-sup test.
Omegatol=sqrt(10000*eps);
% Parameters:
E=1000;
nu=0.24;
uzex= 5.0399e-3;
R= 6;
t=6;
L=t;
ang=35/180*pi;
p=  1;
randshiftmult= 0.37;
% randshiftmult= 0.2;
% parshiftmult= 0.0006113;
parshiftmult= 0.;
% parshiftmult= 0.2;
graphics =   ~true;
scale=200;
nc=1;ny=1;nt=1;
A = [1.44,-0.741,-0.53; -0.626,1.589,-0.913; -0.55,0.43,1.756]+eye(3);
% A = [1.44,-0.741,0.53; -0.626,1.589,-0.913; -0.55,0.43,1.756];
% A = [0.1,-0.2, 0.3; 0.2, 0.1, 0.2; 0.1,-0.3, 0.2]/10;
% A = eye(3);

% mix = 1;
% clear mesd
% mesd(mix).ref=[2:10];
% mesd(mix).ref=[2:6];
% % mesd(mix).ref=[7];
% mesd(mix).ny=1;
% mesd(mix).nc=1;
% mesd(mix).nt=1;
% mix = mix+1;

eix=1;
clear eltyd

eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).ref=[2:6];
eltyd(eix).femmf =@femm_deformation_linear_h8msgso;
eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',2));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
eltyd(eix).styl='k*-';
eix=eix+1;


% eltyd(eix).description='H8R';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).ref=[2:6];
% eltyd(eix).femmf =@femm_deformation_linear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',1));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eltyd(eix).styl='ro-';
% eix=eix+1;


% % 
% eltyd(eix).description='T10';
% eltyd(eix).mf =@T10_block;
% eltyd(eix).ref=[1:4];
% eltyd(eix).femmf =@femm_deformation_linear;
% eltyd(eix).integration_rule=tet_rule (struct('npts',4));
% eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% eltyd(eix).styl='kv--';
% eix=eix+1;
%
% eltyd(eix).description='T10';
% eltyd(eix).mf =@T10_block;
% eltyd(eix).ref=[1:4];
% eltyd(eix).femmf =@femm_deformation_linear;
% eltyd(eix).integration_rule=tet_rule (struct('npts',1));
% eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% eltyd(eix).styl='k^--';
% eix=eix+1;

% eltyd(eix).description='H20R';
% eltyd(eix).mf =@H20_block;
% eltyd(eix).ref=[2:6];
% eltyd(eix).femmf =@femm_deformation_linear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',2));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eltyd(eix).styl='ro-';
% eix=eix+1;

% eltyd(eix).description='NICE-T4';
% eltyd(eix).mf =@T4_block;
% eltyd(eix).femmf =@feblock_defor_ss_nice;
% eltyd(eix).integration_rule= simplex_nq_rule (3);
% eltyd(eix).surface_integration_rule=tri_rule(1);
% eltyd(eix).styl='m^-';
% eix=eix+1;
% %
% eltyd(eix).description='NICE-H8';
% eltyd(eix).mf =@block;
% eltyd(eix).femmf =@feblock_defor_ss_nicer;
% eltyd(eix).integration_rule= tensprod_nq_rule (struct('dim',3,'order',1));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',1));
% eltyd(eix).styl='md-';
% eix=eix+1;

% eltyd(eix).description='NICE-H27';
% eltyd(eix).mf =@blockH27;
% eltyd(eix).femmf =@feblock_defor_ss_nicer;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
% eltyd(eix).styl='mh-';
% eix=eix+1;
%

% eltyd(eix).description='H8R';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@femm_deformation_linear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',1));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eltyd(eix).styl='bx-';
% eix=eix+1;
%
%
% eltyd(eix).description='NICE-H20';
% eltyd(eix).mf =@blockH20;
% eltyd(eix).femmf =@feblock_defor_ss_nicer;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
% eltyd(eix).styl='mo--';
% eix=eix+1;
%
% eltyd(eix).description='H20';
% eltyd(eix).mf =@blockH20;
% eltyd(eix).femmf =@feblock_defor_ss;
% eltyd(eix).integration_rule=gauss_rule (3,3);
% eltyd(eix).surface_integration_rule=gauss_rule(2, 2);
% eltyd(eix).styl='k>--';
% eix=eix+1;

for eix = 1:length(eltyd)
    disp(eltyd(eix).description);
    h=[];
    alpha=[];
    rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
    
    ns=[];
    uzs=[];
    neqns = [];
    
    for     ref=eltyd(eix).ref;
        %% Create the mesh and initialize the geometry
        [fens,fes]= feval (eltyd(eix).mf,t,t,t,nc*ref,ny*ref,nt*ref);
        
        if parshiftmult~=0
            
            ax=parshiftmult*t;
            ay=parshiftmult*t;
            az=parshiftmult*t;
            fens = transform_apply(fens,@(x, data) (x +[-ax/2*x(3)^2*x(2),ay/2*x(3)*x(1)^2,az/4*x(1)^3*x(2)+ax/4*x(3)^2*sin(0.2*x(2)^3)]), []);
        end
        xy=fens.xyz;
        for i=1:count(fens)
            xy(i,:)=xy(i,:)+xy(i,:)*(A);
        end
        fens.xyz=xy;
        bdry_fes = mesh_boundary(fes, struct('other_dimension',1.0));
        
        if (graphics)
            mesh{1}=fens;
            mesh{2}=fes;
            gv=drawmesh(mesh,'fes','facecolor',[0.9,0.95,0.9],'edgecolor','black','shrink', 0.9);
            grid on; font ='Times';
            l=get(gca,'xlabel');
            set (l,'FontSize',14);                set(l,'FontName',font)
            l=get(gca,'ylabel');
            set (l,'FontSize',14);                set(l,'FontName',font)
            l=get(gca,'zlabel');
            set (l,'FontSize',14);                set(l,'FontName',font)
            l=gca;
            set (l,'FontSize',14);                set(l,'FontName',font)
            %                 camset (gv, [])
            pause(1)
        end
        prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
        mater = material_deformation_linear_triax (struct('property',prop ));
        
        
        feb = feval (eltyd(eix).femmf, struct ('material',mater,...
            'fes',fes, 'integration_rule',eltyd(eix).integration_rule));
        % ,'stabfact',0.1));
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        
        %% Define the displacement field, and zero it out
        u   = 0*geom;
        
        %% Apply the EBC''s
        nl=connected_nodes(bdry_fes);
        u = set_ebc(u, nl, nl*0+1, nl*0+1, nl*0);
        u = set_ebc(u, nl, nl*0+1, nl*0+2, nl*0);
        u = set_ebc(u, nl, nl*0+1, nl*0+3, nl*0);
        u = apply_ebc (u);
        
        %% Number the equations
        u   = numberdofs (u);
        
        Gh = infsup_Gh(feb, sysmat_assembler_sparse, geom, u);
        Gh=(Gh+Gh')/2;
        Sh = infsup_Sh(feb, sysmat_assembler_sparse, geom, u);
        Sh=(Sh+Sh')/2;
        [U,S,V] = svd(full(Sh),'econ');
        if (1)
            Ghm =full(Gh);
            clear Gh
            Shm=full (Sh);
            clear Sh
            [modes,Omega]=eig(Ghm,Shm);
        else
            opts.issym= true;
            opts.isreal= true;
            [modes,Omega]=eigs(Gh,Sh,round(0.5*size(Gh,1)),'LM',opts);
        end
        
        % kpm=min([get(u, 'neqns'),(3*get(u, 'neqns')-get(u, 'neqns')+1)+2 ])
        %             kpm =round(get(u, 'neqns')/8)
        %             [modes,Omega]=eigs(get(Gh,'mat'),get(Sh,'mat'),kpm,'SM');
        clear Ghm Shm
        absOmega =diag(Omega);
        ix=find(absOmega<0);
        absOmega(ix)=0;
        absOmega =(sort(real(absOmega)))';
        %                     loglog(absOmega,eltyd(eix).styl); hold on; pause(1);
        ix=find(absOmega>Omegatol);
        if isempty(ix)
            error (' Bad guess of the number of eigenvalues')
        end
        %         rank(Shm)-rank(Ghm)
        %         ix(1)
        format short
        alpha=[alpha,absOmega(ix(1))]
        h=[h,1/(count(fens))^(1/3)];
        save([eltyd(eix).description,'-ref',num2str(ref),'-results.mat'],'absOmega','h','alpha')
    end
    loglog(h,alpha,eltyd(eix).styl); hold on
    %         set(gca,'YLim', [0.0001, 0.1]);
    font='Times';
    xlabel ('1/N')
    ylabel ('\alpha')
    h=legend(eltyd(eix).description,'Location','Southeast');
    set (h,'FontSize', 16);
    set(gca,'Position',[0.15 0.15 0.5 0.5]);
    set(gca,'FontSize',14);
    set(gca,'FontName',font)
    l=get(gca,'xlabel');
    set (l,'FontSize',16);
    set(l,'FontName',font)
    l=get(gca,'ylabel');
    set (l,'FontSize',16);
    set(l,'FontName',font)
    hold on; grid on; figure (gcf);
    %         print(gcf, '-r0', [mfilename '-' eltyd(eix).description'.eps'],'-deps2c');
    figure(gcf)
    pause(1)
%         close (gcf)
end
