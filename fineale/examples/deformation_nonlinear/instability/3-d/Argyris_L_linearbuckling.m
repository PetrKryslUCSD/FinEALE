%
function Argyris_L_linearbuckling
E = 71240;
nu=.31;
%     nu=.0;
stabfact = 0*E;
rho=5e-9;
graphics = ~false;
Cam =1.0e+003 * [  -1.4407   -0.8522    1.1901    0.1089    0.1350   -0.0492         0         0    0.0010    0.0062];
scale=1;
nw = 2;
nL = 1;
nt = 1;
L=255;
w=30;
t=0.6;
Fmag= 1e-5;

rand('state',0);% try to comment out this line and compare
%                   results for several subsequent runs
%
    function [fens,fes] =  Lframe_Mesh(L,w,t,nw,nL,nt)
        Nodes = [
            1    L-w,0,0
            2    L,0,0
            3     L, L,0
            4     L-w,L-w,0
            5     0,L,0
            6     0,L-w,0
            ];
        [fens1,fes1] = Q4_quadrilateral(Nodes( [1, 2, 3, 4],2:4),nw,nL,struct('other_dimension',1));
        [fens2,fes2] = Q4_quadrilateral(Nodes( [4,3,5,6],2:4),nw,nL,struct('other_dimension',1));
        [fens,g1,g2] = merge_meshes(fens1, fes1, fens2, fes2, 0.01);
        fes= cat(g1,g2);
        [fens,fes] = H8_extrude_Q4(fens,fes,nt,@(xyz, layer)([xyz(1:2), (t/nt)*(layer-1)]));
    end
    function [fens,fes] = h8H8(nw,nL)
        [fens,fes] = Lframe_Mesh(L,w,t,nw,nL,2*nt);
    end
    function [fens,fes] = h8H27(nw,nL)
        [fens,fes] = Lframe_Mesh(L,w,t,nw,nL,nt);
        [fens,fes] = H8_to_H27(fens,fes);
    end
    function [fens,fes] = h8H20(nw,nL)
        [fens,fes] = Lframe_Mesh(L,w,t,nw,nL,nt);
        [fens,fes] = H8_to_H20(fens,fes);
    end
    function [fens,fes] = h8H64(nw,nL)
        [fens,fes] = Lframe_Mesh(L,w,t,nw,nL,nt);
        [fens,fes] = H8_to_H64(fens,fes);
    end

eix=1;
clear eltyd



eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@h8H8;
eltyd(eix).blf =@femm_deformation_nonlinear_h8msgso;
eltyd(eix).integration_rule=gauss_rule(struct('dim',3,'order',2));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
eltyd(eix).styl='ko-';
eix=eix+1;

% eltyd(eix).description='H8';
% eltyd(eix).mf =@h8H8;
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
% eltyd(eix).mf =@h8H20;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3, 'order',2));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
% eltyd(eix).styl='ks-';
% eix=eix+1;

    function [Bucklen, Bucklep] =Compute
        %          Mesh
        [fens,fes] = eltyd(eix).mf (nw, nL);
        %         [fens, fes]=renum_mesh(fens, fesAnd, 'symrcm');
        % Material
        prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
        mater = material_deformation_stvk_triax (struct('property',prop ));
        femm = eltyd(eix).blf (struct ('material',mater, 'fes',fes,'integration_rule',eltyd(eix).integration_rule));
        bfes =mesh_boundary(fes, []);
        icl = fe_select(fens, bfes, struct('box', [-Inf,Inf,0,0,-Inf,Inf],'inflate',1/100));
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
        %         Bottom surface
        ebc_fenids=fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',1/1000));
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
        Traction=zeros(3,1); Traction(1) =Fmag/t/w;
        fi=force_intensity(struct('magn', Traction));
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
            draw(femm,gv, struct ('x', geom,'u',+100/abs(Fmag)*u));
            draw(femm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
        end
        
        Kgeo = stiffness_geo(femm1, sysmat_assembler_sparse, geom, u,0*u);
        
        neigvs = 2;
        options.tol =1e-206;
        options.maxit = 5000;
        options.disp = 0;
        [W,Omega]=eigs(-(Kgeo+Kgeo')/2,(K+K')/2,...
            neigvs,'lm', options);
        %         diag(Omega)
        %         [Omegas,ix]=sort(1./diag(Omega)) ;
        Omegas=(1./(diag(Omega)));
        Bucklen =Omegas(1); Bucklep=Omegas(2);
        % Plot
        scale=100;
        for i=(1:neigvs)
            phi = scatter_sysvec(u,W(:,i));
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
                gv=reset (gv,[]);
                %                 set(gca,'FontSize',16)
                %     camset(gv,[-0.1591    0.0085   -0.3252    0.0001    0.0043    0.0013   -0.8985    0.0237    0.4383    5.4322]);
                draw(femm,gv, struct ('x', geom,'u',0*phi, 'facecolor','none'));
                draw(femm,gv, struct ('x', geom,'u',+scale*phi,'colorfield',colorfield));
                camset(gv,1.0e+003 *[1.9128   -1.0031    0.6488    0.1593    0.1575   -0.0344         0         0    0.0010    0.0048]);
                %     text(0,1.1*R,1.1*R,['\omega_' num2str(i) '=' num2str(sqrt(Omega(i,i))/2/pi) ],'FontSize',24);
                axis off; set_graphics_defaults
                %                 pause
            end
        end
    end

for eix = 1:length(eltyd)
    Bucklens=[]; Buckleps = [];
    nelperedge =(1:9);
    for nL =  nelperedge
        [Bucklen, Bucklep] = Compute;
        Bucklens=[Bucklens,Bucklen];
        Buckleps = [Buckleps,Bucklep];
    end
    disp(['Data{end+1}=[']);
    disp([num2str(nelperedge,15)]);
    disp([num2str(Bucklens,15)]);
    disp([num2str(Buckleps,15)]);
    disp(['];' 'Style{end+1}=''' name_to_style(eltyd(eix).description) ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
    %         load_displacement
    %         Bucklens, Buckleps
end
end