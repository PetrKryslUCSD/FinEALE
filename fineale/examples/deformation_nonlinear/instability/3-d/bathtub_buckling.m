function bathtub_buckling
E = 448.88;
%     [x,y]=solve ('400889.8-x/3/( 1-2*y)','80.194-x/2/(1+y)')
nu= 0.49;
nu= 0.4999999;
%     nu= 0.3;
%    nu=.0;
rho=5e-9;
graphics = ~false;
Cam =1.0e+003 * [-1.4407   -0.8522    1.1901    0.1089    0.1350   -0.0492         0         0    0.0010    0.0062];
scale=1;
mult =3;
nw =mult*3;
nL = mult*1;
%     nw =40;
%     nL = 20;
nt = 1;
L=16;
w= 52;
t=1;
Fmag= 1/t/w;
Pcritt =pi^2*E*(t^3*w/12)/(2*L)^2;
Pcritw =pi^2*E*(t*w^3/12)/(2*L)^2;

prop = property_deformation_neohookean (struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));

%     prop = property_deformation_linear_iso(struct('E',E,'nu',nu));
%     mater = material_deformation_linear_triax(struct('property',prop));
%

rand('state',0);% try to comment out this line and compare
%                   results for several subsequent runs
%
    function [fens,fes] = h8T4(nw,nL)
    [fens,fes]=t4block(w,t,L, nw,nt,nL);
    end
    function [fens,fes] = h8T4d(nt)
    [fens,fes]=t4blocka(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
    [fens1,fes1]=t4blockb(1.0, 1.0, 1.0, nt(1),nt(2),nt(3));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, 1/1000);
    fes=cat(2,fes1,fes2);
    end
    function [fens,fes] = h8H8(nw,nL)
    [fens,fes]=H8_block(w,t,L, nw,nt,nL);
    end
    function [fens,fes] = h8H27(nw,nL)
    [fens,fes]=H27_block(w,t,L, nw,nt,nL);
    end
    function [fens,fes] = h8H20(nw,nL)
    [fens,fes]=block(w,t,L, nw,nt,nL);
    [fens,fes] = H8_to_H20(fens,fes);
    end
    function [fens,fes] = h8H64(nw,nL)
    [fens,fes]=block(w,t,L, nw,nt,nL);
    [fens,fes] = H8_to_H64(fens,fes);
    end

eix=1;
clear eltyd

eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@h8H8;
eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgso(struct ('material',mater, ...
    'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
eix=eix+1;

% eltyd(eix).label ='T4';
% eltyd(eix).mf =@h8T4;
% eltyd(eix).blf =@femmlock_defor_ss;
% eltyd(eix).integration_rule= tet_rule(1);
% eltyd(eix).surface_integration_rule=tri_rule(1);
% eltyd(eix).styl='rv-.';
% eix=eix+1;
%
% eltyd(eix).label ='H8';
% eltyd(eix).mf =@h8H8;
% eltyd(eix).blf =@femmlock_defor_nonlinear;
% eltyd(eix).integration_rule= gauss_rule (3,2);
% eltyd(eix).surface_integration_rule=gauss_rule (2,2);
% eltyd(eix).styl='rx-.';
% eix=eix+1;

%         eltyd(eix).label ='NICE-H8';
%         eltyd(eix).mf =@h8H8;
%         eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
%         eltyd(eix).integration_rule= tensprod_nq_rule (struct('dim',3,'order',1));
%         eltyd(eix).surface_integration_rule=tensprod_nq_rule(struct('dim',2,'order',1));
%         eltyd(eix).styl='rx-.';
%         eix=eix+1;

%         eltyd(eix).label ='NICE-T4';
%         eltyd(eix).mf =@h8T4;
%         eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
%         eltyd(eix).integration_rule= simplex_nq_rule (3);
%         eltyd(eix).surface_integration_rule=tri_rule(1);
%         eltyd(eix).styl='m^-';
%         eix=eix+1;


%     eltyd(eix).label ='H27';
%     eltyd(eix).mf =@h8H27;
%     eltyd(eix).blf =@(fes)femm_deformation_nonlinear(struct ('material',mater, ...
%         'fes',fes, ...
%         'integration_rule',gauss_rule(struct('dim',3,'order',3))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',3));
%     eltyd(eix).styl='ks-';
%     eix=eix+1;

%     eltyd(eix).label ='H64';
%     eltyd(eix).mf =@h8H64;
%     eltyd(eix).blf =@femmlock_defor_nonlinear;
%     eltyd(eix).integration_rule=gauss_rule (3,4);
%     eltyd(eix).surface_integration_rule=gauss_rule(2, 4);
%     eltyd(eix).styl='ks-';
%     eix=eix+1;

%         eltyd(eix).label ='NICE-H27';
%         eltyd(eix).mf =@h8H27;
%         eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
%         eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
%         eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
%         eltyd(eix).styl='mh-';
%         eix=eix+1;

%     eltyd(eix).label ='H20R';
%     eltyd(eix).mf =@h8H20;
%     eltyd(eix).blf =@femmlock_defor_nonlinear;
%     eltyd(eix).integration_rule=gauss_rule (3,3);
%     eltyd(eix).surface_integration_rule=gauss_rule(2, 2);
%     eltyd(eix).styl='ro-';
%     eix=eix+1;
%

% eltyd(eix).label ='NICE-H64';
% eltyd(eix).mf =@h8H64;
% eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',3));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',3));
% eltyd(eix).styl='mh-';
% eix=eix+1;

    function Compute
    %          Mesh
    [fens,fes] = eltyd(eix).mf (nw, nL);
    femm = eltyd(eix).blf (fes);
    bfes =mesh_boundary(fes, []);
    icl = fe_select(fens, bfes, struct('box', [-Inf,Inf,-Inf,Inf,L,L],'inflate',1/100));
    %          gv=drawmesh( {fens,fes},'facecolor','none'); gv=drawmesh( {fens,bfes(icl)},'gv',gv,'facecolor','blue');
    %         view(3); pause(1)
    
    efemm = femm_deformation(struct ('material',mater, 'fes',subset(bfes,icl),...
        'integration_rule',eltyd(eix).surface_integration_rule));
    sfemm = femm_deformation(struct ('material',mater, 'fes',bfes,...
        'integration_rule',eltyd(eix).surface_integration_rule));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = 0*geom; % zero out
    % Apply EBC's
    %         Bottom surface
    ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0],'inflate',1/1000));
    ebc_prescribed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fenids*0+3;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    ebc_fenids=fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',1/1000));
    ebc_prescribed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fenids*0+1;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    ebc_fenids=fenode_select (fens,struct('box',[w,w,-Inf,Inf,-Inf,Inf],'inflate',1/1000));
    ebc_prescribed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fenids*0+1;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    ebc_fenids =(1:count(fens));
    ebc_prescribed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fenids*0+2;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    
    u   = apply_ebc (u);
    
    % Number equations
    u   = numberdofs (u);
    u = u*0; % zero out the displacement
    
    
    % Update the FEMM
    femm  =update(femm,geom,u,u);
    femm1=femm;                                           % Make a copy of the state
    Load=zeros(3,1); Load(3) =-Fmag;
    fi=force_intensity(struct('magn',Load));
    F = distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
    %         sum (get(F,'vec'))
    
    K = stiffness(femm1, sysmat_assembler_sparse, geom, u,u);
    u = scatter_sysvec(u, K\F); % Displacement increment
    %         gv=graphic_viewer;
    %         gv=reset (gv,[]);
    %         draw(femm,gv, struct ('x', geom,'u',+1e3*u));
    %         draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
    
    Kgeo = stiffness_geo(femm1, sysmat_assembler_sparse, geom, u, 0*u);
    
    neigvs = 2;
    options.tol =1e-206;
    options.maxit = 5000;
    options.disp = 0;
    [W,Omega]=eigs(-(Kgeo+Kgeo')/2,(K+K')/2,neigvs,'lm', options);
    diag(Omega)
    %         [Omegas,ix]=sort(1./diag(Omega)) ;
    Omegas=(1./(diag(Omega)))
    
    Control=fenode_select (fens,struct('box',[0, 0,0,0,L,L],'inflate',1/1000));
    Controleq =gather_dofnums(u, Control);
    
    for i=(1:neigvs)
        phi = scatter_sysvec(u,sign(W(Controleq(3),i))*W(:,i));
        phiv = phi.values;
        norm(K*gather_sysvec(phi) +Omegas(i)*Kgeo*gather_sysvec(phi))
        disp(['  Buckling load '  num2str((Omegas(i)*Fmag)) ]);
        %     clf;
        if (graphics)
            gv=graphic_viewer;
            gv=reset (gv,[]);
            set_graphics_defaults
            %                 set_label_defaults
            scale =L/6/max(abs(gather_sysvec(phi)));
            phivmag = sqrt(phiv(:,1).^2+phiv(:,2).^2+phiv(:,3).^2);
            dcm=data_colormap(struct ('range',[min(phivmag),max(phivmag)], 'colormap',cadcolors));
            colors=map_data(dcm, phivmag);
            colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
            gv=reset (gv,[]);
            %                 set(gca,'FontSize',16)
            %     camset(gv,[-0.1591    0.0085   -0.3252    0.0001    0.0043    0.0013   -0.8985    0.0237    0.4383    5.4322]);
            %                 draw(femm,gv, struct ('x', geom,'u',0*phi, 'facecolor','none'));
            %                 draw(femm,gv, struct ('x', geom,'u',+scale*phi,'colorfield',colorfield));
            draw(femm,gv, struct ('x', geom,'u',scale*phi, 'facecolor','white'));
            %                 set(gcf,'Position', [261   219   431   213]);
            camset(gv,[26.7898 -449.8332   18.8926   26.7898    0.5000   18.8926         0         0    1.0000 4.1547 ]);
            axis off
            text(0,1.1*w,1.3*L,['\lambda_{cr}= ' num2str(Omegas(i)*Fmag) ],'FontSize',34);
            saveas (gcf,[mfilename '-' eltyd(eix).description '-nw' num2str(nw) '-m' num2str(i)  '.eps']);
            %                printeps (gcf,[mfilename '-' eltyd(eix).label '-nw' num2str(nw) '-m' num2str(i) '-s' num2str(stabfact) '.png']);
        end
    end
    end

for eix = 1:length(eltyd)
    %         for stabfact= [0, 0.01,0.05, 0.1]
    Compute;
end
end