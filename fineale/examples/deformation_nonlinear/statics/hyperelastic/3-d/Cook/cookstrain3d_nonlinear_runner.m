function cookstrain3dnonlinear_runner
graphics= ~true;
doeig=~true;
aspect=1;
smult=0;
ne= [1,2];
%     nu =1/2*lambda/(G+lambda)
nu= 0.4999;
%     nu= 0.49999;
%     nu= 0.3;
%     E=G*(3*lambda+2*G)/(G+lambda);
E=240.565;%MPa
magn =100/16;
nincr =29;
% Data: Elguedj et al. CMAME 2008
%
convutip=6.9083;% for 128x128 mesh

prop = property_deformation_neohookean (struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));

clear eltyd
eix=1;
 
eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).convertf =[];
eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
eltyd(eix).ne=16:-2:2;
%               eltyd(eix).ne= 1;
eix=eix+1;

% eltyd(eix).description='H20R';
% eltyd(eix).mf =@H8_block;
% eltyd(eix).convertf =@H8_to_H20;
% eltyd(eix).blf =@(fes)femm_deformation_nonlinear(struct ('material',mater, 'fes',fes, ...
%     'integration_rule',gauss_rule(struct('dim',3,'order',2))));;
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',4));
% eltyd(eix).ne= 2:2:9;
% eix=eix+1;

%    eltyd(eix).description='H27';
%    eltyd(eix).mf =@H8_block;
%    eltyd(eix).convertf =@H8_to_H27;
%    eltyd(eix).blf =@(fes)femm_deformation_nonlinear(struct ('material',mater, 'fes',fes, ...
%                     'integration_rule',gauss_rule(struct('dim',3,'order',3))));;
%    eltyd(eix).integration_rule= gauss_rule(struct('dim',3,'order',3));
%    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',4));
%    eltyd(eix).ne= 8:-1:1;
%    eix=eix+1;

%    eltyd(eix).description='H64';
%    eltyd(eix).mf =@H8_block;
%    eltyd(eix).convertf =@H8_to_H64;
%    eltyd(eix).blf =@femm_deformation_nonlinear;
%    eltyd(eix).integration_rule= gauss_rule(struct('dim',3,'order',4));
%    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',4));
%    eltyd(eix).ne= 8%1:-1:1;
%    eix=eix+1;
%
%    eltyd(eix).description='H20';
%    eltyd(eix).mf =@H8_block;
%    eltyd(eix).convertf =@H8_to_H20;
%    eltyd(eix).blf =@femm_deformation_nonlinear;
%    eltyd(eix).integration_rule= gauss_rule(struct('dim',3,'order',2));
%    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
%    eltyd(eix).ne= 8%1:-1:1;
%    eix=eix+1;



for eix = 1:length(eltyd)
    
    % Mesh
    h=[];
    utip=[];
    
    %     for n=[4,8,16,32]
    for n=eltyd(eix).ne
        % for n=[2,4,8,16,32,64,128]
        %             for n=[2,4,8,16,24,32,64,128]
        %         for n=[4]
        options =struct('other_dimension', 1.0);
        [fens,fes] = eltyd(eix).mf(48,44,16/n, n, aspect*n, 1);
        dxy=min(48,44)/n/aspect/100;
        sxy=smult*dxy;
        rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
        xys=fens.xyz;
        for i=1:count(fens)
            xy=xys(i,:);
            if (xy(1)>dxy) & (xy(2)>dxy) & (xy(1)<48-dxy) & (xy(2)<44-dxy)
                xy=[xy(1)+(2*rand(1)-1)/2*sxy,xy(2)+(2*rand(1)-1)/2*sxy,xy(3)];
            end
            xys(i,:)=xy;
        end
        fens.xyz=xys;
        xys=fens.xyz;
        for i=1:count(fens)
            xy=xys(i,:);
            y=xy(2) + (xy(1)/48)*(44 -xy(2)/44*28);
            xys(i,:)=[xy(1) y xy(3)];
        end
        fens.xyz=xys;
        if (~isempty( eltyd(eix).convertf ))
            [fens,fes] = eltyd(eix).convertf(fens,fes);
        end
        %             mesh{1}=fens;
        %             mesh{2}=gcells;
        %             drawmesh(mesh); view(2); pause (1)
        femm = eltyd(eix).blf(fes);
        
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = clone(geom,'u');
        u   = u*0; % zero out
        % Apply EBC's
        maxz =max(fens.xyz(:,3));
        ebc_fenids=fenode_select (fens,struct('box',[0,0,0, 44, 0,maxz],'inflate',44/10000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=[];
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        ebc_fenids=(1:count(fens));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+3;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        % Number equations
        u   = numberdofs (u);
        % Boundary
        bdry_fes = mesh_boundary(fes, struct('other_dimension', 1.0));
        icl = fe_select(fens, bdry_fes, struct('box', [48,48,44,60,0,maxz],'inflate',48/1000)) ;
        efemm = femm_deformation(struct ('material',mater,'fes', subset(bdry_fes,icl),...
            'integration_rule',eltyd(eix).surface_integration_rule));
        % Solve
        % Now comes the nonlinear solution
        tup = 1;
        u = u*0; % zero out the displacement
        utol = 1e-10*u.nfreedofs;
        
        if (graphics),gv=graphic_viewer;end
        
        t=0; % time=load magnitude
        incr=1;
        % Update the FEMM
        femm  =update(femm,geom,u,u);
        while (incr <= nincr)
            t = t + tup / nincr;
            disp(['Increment ' num2str(incr) ]); % pause
            % Initialization
            u1 = u; % guess
            u1 = apply_ebc(u1);
            du = 0*u; % this will hold displacement increment
            du = apply_ebc(du);
            iter=1;
            %             gv=reset(clear(gv,[]),[]);
            %             draw(femm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
            femm1=femm;                                           % Make a copy of the state
            while 1 % Iteration loop
                fi=force_intensity(struct('magn',[0;t*magn;0]));
                F = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);% cconfiguration-independent loads
                F = F + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
                K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
                % Displacement increment
                du = scatter_sysvec(du, K\F);
                u1 = u1 + du;                                       % increment displacement
                %                 draw(femm,gv, struct ('x', geom,'u',u1, 'facecolor','none'));
                %                 pause(1)
                disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
                if (norm(du) < utol) break; end;                    % convergence check
                iter=iter+1;
            end %while 1 % Iteration loop
            % Update the FEMM
            femm1  =update(femm1,geom,u1,u);
            u = u1;                                               % update the displacement
            if (graphics)
                gv=reset(clear(gv,[]),[]);
                cmap = jet(16);
                %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
                %             draw(sfemm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
                %                     fld = field_from_integration_points(femm, geom, u1, [], 'pressure',1);
                %                     nvals=get(fld,'values');%min(nvals),max(nvals)
                %                     nvalsrange=[min(nvals),max(nvals)];
                %                     dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
                %                     colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
                %                     draw(femm,gv, struct ('x', geom, 'u', u,'colorfield',colorfield, 'shrink',1.0));
                draw(femm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
                draw(femm,gv, struct ('x', geom, 'u', u1,'facecolor','yellow', 'shrink',1.0));
                % draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
                %                     colormap(cmap);
                %                     cbh=colorbar;
                %                     set(cbh,...
                %                         'Position',[0.815 0.15 0.05 0.7],...
                %                         'YLim',[0,1],...
                %                         'YTick',[0,1],...
                %                         'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
                %                     set(get(cbh,'XLabel'),'String','\sigma_z');
                pause(0.5);
            end
            incr = incr + 1;
        end %while (incr <= nincr)
        
        % Extract the solution
        nl=fenode_select (fens,struct ('box',[48,48,60,60,0,1],'inflate',2/100)) ;
        theutip=gather_values(u,nl);
        theutip =sum(theutip(:,2))/length(nl);
        disp ([' el/edge=' num2str(n) ', displ=' num2str(theutip, 15)])
        h=[h 44/n];
        utip =[utip theutip];
        
        % Plot
        if graphics
            gv=graphic_viewer;
            gv=reset (gv,struct ([]));
            scale=1;
            cmap = gray;
            cmap=jet;
            %             fld = field_from_integration_points(femm, geom, u, [], 'vol_strain', 1);
            %             fld = field_from_integration_points(femm, geom, u, [], 'Cauchy', 1);
            fld = field_from_integration_points(femm, geom, u, [], 'pressure', 1);
            nvals=fld.values;%min(nvals),max(nvals)
            %             nvalsrange=[-0.05, -0.02];%[min(nvals),max(nvals)]
            nvalsrange=[min(nvals),max(nvals)];
            dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
            draw(femm,gv, struct ('x', geom, 'u', scale*u,'colorfield',colorfield,...
                'edgecolor','white','shrink',1.0));
            xlabel('X','FontName','Times')
            ylabel('Y','FontName','Times')
            set(gca,'FontSize',14);
            set(gca,'FontName','Times')
            % geom3d=combine (geom, 0*fld);
            % u3d=combine (u, 100000000*fld);
            % draw(femm,gv, struct ('x', geom3d, 'u', +scale*u3d,'colorfield',colorfield, 'shrink',1.0));
            colormap(cmap);
            cbh=colorbar;
            set(cbh,'FontSize',14);
            set(cbh,'FontName','Times')
            set(cbh,...
                'Position',[0.625 0.247 0.035 0.51],...
                'YLim',[0,1],...
                'YTick',[0,1],...
                'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
            %             set(get(cbh,'XLabel'),'String','Volumetric strain');
            view (2)
            
        end
    end
    disp(eltyd(eix).description);
    disp([h;utip]);
    %     [xestim, beta, residual] =
    %     richextrapol(utip(end-2:end),h(end-2:end))
end
end
