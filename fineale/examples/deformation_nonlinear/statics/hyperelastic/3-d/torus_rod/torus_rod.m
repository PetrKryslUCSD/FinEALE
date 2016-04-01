function   torus_rod
E=3*5.67e6;
nu=0.499;
%     E=G*(3*lambda+2*G)/(G+lambda);
rho=5e-9;
R= 1;
sR =9;
rA=pi/8;
t= pi/2;
Fmag=0.195e6;

nincr = 5;
graphics = ~false;
cmap = jet;
scale=1;
nR = [2];
nt = 2*[4];

%
    function [fens] = to_cyl_surf(fens,fes)
        bdry_fes = mesh_boundary(fes, []);
        tol=0.00001*t;
        bcl = [fe_select(fens, bdry_fes, ...
            struct ('box',[0 0 -100*R 100*R -100*R 100*R],'inflate',tol/2)),...
            fe_select(fens, bdry_fes, ...
            struct ('box',[t t -100*R 100*R -100*R 100*R],'inflate',tol/2)),...
            fe_select(fens, bdry_fes, ...
            struct ('box',[0 t 0 0 -100*R 100*R],'inflate',tol/2))];
        conns=bdry_fes.conn;
        xyzs=fens.xyz;
        for j=setdiff((1:count(bdry_fes)),bcl)
            conn=conns(j,:);
            for k=1:length(conn)
                xyz=xyzs(conn(k),:);
                xyzs(conn(k),:)=[xyz(1),R/norm(xyz(2:3))*xyz(2:3)];
            end
        end
        fens.xyz=xyzs;
    end
    function [fens,fes] = t4cyl(t,R, nt,nR)
        [fens,fes] = h8cylH8(t,R, nt,nR);
        [fens,fes] = t4del(fens);
    end
    function [fens,fes] = t4cylT10(t,R, nt,nR)
        [fens,fes] = T4_cylinderdel(t,R, nt,nR);
        [fens,fes] = T4_to_T10(fens,fes);
        [fens] = to_cyl_surf(fens,fes);
    end
    function [fens,fes] = h8cylH8(t,R, nt,nR)
        [fens,fes]=Q4_circle(R,log2(nR),1.0);
        [fens1,fes1] = mirror_mesh(fens, fes, [1,0], [0,0]);
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, t/10000);
        fes=cat(fes1,fes2);
        [fens,fes] = H8_extrude_Q4(fens,fes,nt,@(xyz,k)([xyz,k*t/nt]));
        fens = transform_apply(fens,@(x,d)([x(3),x(2),-x(1)]), []);
    end
    function [fens,fes] = scylH8(t,R, nt,nR)
        fens(1) =fenode(struct('id',1,'xyz',[0,0]));
        fens(2) =fenode(struct('id',2,'xyz',[0,R]));
        fens(3) =fenode(struct('id',3,'xyz',[R*cos(3*pi/8),R*sin(3*pi/8)]));
        fens(4) =fenode(struct('id',4,'xyz',[R*cos(pi/8),R*sin(pi/8)]));
        fens(5) =fenode(struct('id',5,'xyz',[R,0]));
        fens(6) =fenode(struct('id',6,'xyz',[R/2,0]));
        %         [fens,fes]=q4circle(R,log2(nR),1.0);
        fes(1)=fe_Q4(struct('id',1,'conn',[1,6,3,2]));
        fes(2)=fe_Q4(struct('id',2,'conn',[6,5,4,3]));
        [fens1,fes1] = mirror_mesh(fens, fes, [1,0], [0,0]);
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, t/10000);
        fes=cat(2,fes1,fes2);
        %         [fens1,fes1] = mirror_mesh(fens, fes, [0,1], [0,0]);
        %         [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, t/10000);
        %         fes=cat(2,fes1,fes2);
        [fens,fes] = extrudeq4(fens,fes,nt,@(xyz,k)([xyz,k*t/nt]));
        fens = transform_apply(fens,@(x,d)([x(3),x(2),-x(1)]), []);
    end
    function [fens,fes] = s1cylH8(t,R, nt,nR)
        [fens,fes]=q4circle(R,log2(nR),1.0);
        fens(4) =fenode(struct('id',4,'xyz',[R*cos(pi/8),R*sin(pi/8)]));
        fens(7) =fenode(struct('id',7,'xyz',[R/2*cos(pi/8),R/2*sin(pi/8)]));
        [fens1,fes1] = mirror_mesh(fens, fes, [1,0], [0,0]);
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, t/10000);
        fes=cat(2,fes1,fes2);
        %         [fens1,fes1] = mirror_mesh(fens, fes, [0,1], [0,0]);
        %         [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, t/10000);
        %         fes=cat(2,fes1,fes2);
        [fens,fes] = extrudeq4(fens,fes,nt,@(xyz,k)([xyz,k*t/nt]));
        fens = transform_apply(fens,@(x,d)([x(3),x(2),-x(1)]), []);
    end
    function [fens,fes] = h8cylH27(t,R, nt,nR)
        [fens,fes] = h8cylH8(t,R, nt,nR);
        [fens,fes] = H8_to_H27(fens,fes);
        [fens] = to_cyl_surf(fens,fes);
    end
    function [fens,fes] = scylH27(t,R, nt,nR)
        [fens,fes] = scylH8(t,R, nt,nR);
        [fens,fes] = H8_to_H27(fens,fes);
        [fens] = to_cyl_surf(fens,fes);
    end
    function [fens,fes] = scylH64(t,R, nt,nR)
        [fens,fes] = scylH8(t,R, nt,nR);
        [fens,fes] = H8_to_H64(fens,fes);
        [fens] = to_cyl_surf(fens,fes);
    end
    function [fens,fes] = s1cylH64(t,R, nt,nR)
        [fens,fes] = s1cylH8(t,R, nt,nR);
        [fens,fes] = H8_to_H64(fens,fes);
        [fens] = to_cyl_surf(fens,fes);
    end

eix=1;
clear eltyd
    
eltyd(eix).description='H8MSGS';
eltyd(eix).mf =@h8cylH8;
eltyd(eix).blf =@femm_deformation_nonlinear_h8msgs;
eltyd(eix).integration_rule= gauss_rule(struct('dim',3,'order',2));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',4));
eix=eix+1;
     
% eltyd(eix).label ='T4';
% eltyd(eix).mf =@t4cyl;
% eltyd(eix).blf =@femmlock_defor_ss;
% eltyd(eix).integration_rule= tet_rule(1);
% eltyd(eix).surface_integration_rule=tri_rule(1);
% eltyd(eix).styl='rv-.';
% eix=eix+1;
%
% eltyd(eix).label ='H8';
% eltyd(eix).mf =@h8cylH8;
% eltyd(eix).blf =@femmlock_defor_nonlinear;
% eltyd(eix).integration_rule= gauss_rule (3,2);
% eltyd(eix).surface_integration_rule=gauss_rule (2,2);
% eltyd(eix).styl='rx-.';
% eix=eix+1;
%
% eltyd(eix).label ='NICE-T4';
% eltyd(eix).mf =@t4cyl;
% eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
% eltyd(eix).integration_rule= simplex_nq_rule (3);
% eltyd(eix).surface_integration_rule=tri_rule(1);
% eltyd(eix).styl='m^-';
% eix=eix+1;
%
% eltyd(eix).label ='NICE-H8';
% eltyd(eix).mf =@h8cylH8;
% eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
% eltyd(eix).integration_rule= tensprod_nq_rule (struct('dim',3,'order',1));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',1));
% eltyd(eix).styl='md-';
% eix=eix+1;

% eltyd(eix).label ='T10';
% eltyd(eix).mf =@t4cylT10;
% eltyd(eix).blf =@femmlock_defor_nonlinear;
% eltyd(eix).integration_rule=tet_rule (4);
% eltyd(eix).surface_integration_rule=tri_rule(3);
% eltyd(eix).styl='kv--';
% eix=eix+1;
%
% eltyd(eix).label ='H27';
% eltyd(eix).mf =@h8cylH27;
% eltyd(eix).blf =@femm_deformation_nonlinear;
% eltyd(eix).integration_rule=gauss_rule(struct('dim',3,'order',3));;
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',4));;
% eltyd(eix).styl='ks-';
% eix=eix+1;

% eltyd(eix).label ='NICE-H27';
% eltyd(eix).mf =@scylH27;
% eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
% eltyd(eix).styl='mh-';
% eix=eix+1;

% eltyd(eix).label ='NICE-H64';
% eltyd(eix).mf =@s1cylH64;
% eltyd(eix).blf =@femmlock_defor_nonlinear_nice;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',3));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',3));
% eltyd(eix).styl='mh-';
% eix=eix+1;

% eltyd(eix).label ='H20R';
% eltyd(eix).mf =@uh8cylH20;
% eltyd(eix).blf =@femmlock_defor_ss;
% eltyd(eix).integration_rule=gauss_rule (3,2);
% eltyd(eix).surface_integration_rule=gauss_rule(2, 2);
% eltyd(eix).styl='ro-';
% eix=eix+1;
% %
%
%
% eltyd(eix).label ='NICE-H20';
% eltyd(eix).mf =@uh8cylH20;
% eltyd(eix).blf =@femmlock_defor_ss_nice;
% eltyd(eix).integration_rule=tensprod_nq_rule (struct('dim',3,'order',2));
% eltyd(eix).surface_integration_rule=tensprod_nq_rule (struct('dim',2,'order',2));
% eltyd(eix).styl='mo--';
% eix=eix+1;
% %
% eltyd(eix).label ='H20';
% eltyd(eix).mf =@uh8cylH20;
% eltyd(eix).blf =@femmlock_defor_ss;
% eltyd(eix).integration_rule=gauss_rule (3,3);
% eltyd(eix).surface_integration_rule=gauss_rule(2, 2);
% eltyd(eix).styl='k>--';
% eix=eix+1;

    function [neqns,energy] =Compute(nR,nt)
        nus = [(1-1.8.^(-(1:14)))*nu,(1:100)*0+nu];
        %          Mesh
        [fens,fes] = eltyd(eix).mf (t,R, nt,nR);
        bfes =mesh_boundary(fes, []);
        icl = fe_select(fens, bfes, struct('box', [0,rA,-R*sin(rA),R*sin(rA),R*cos(rA),Inf],'inflate',R/100));
        xyzs=fens.xyz;
        for i=1:size(xyzs,1)
            xyz=xyzs(i,:);
            a=xyz(1); y=xyz(2); z=xyz(3);
            xyz= [0, y, z+sR];
            xyz = xyz*rotmat([0, -a, 0]);
            xyzs(i,:)=xyz;
        end
        fens.xyz=xyzs;
        %         gv=drawmesh( {fens,fes},'facecolor','none'); gv=drawmesh( {fens, subset(bfes,icl)},'gv',gv,'facecolor','blue');
        %         view(3); pause(1)
        % Material
        prop = property_deformation_neohookean (struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));
        % Finite element block
        femm = eltyd(eix).blf (struct ('material',mater, 'fes',fes,...
            'integration_rule',eltyd(eix).integration_rule));
        efemm = eltyd(eix).blf (struct ('material',mater, 'fes',  subset(bfes,icl),...
            'integration_rule',eltyd(eix).surface_integration_rule));
        sfemm = eltyd(eix).blf (struct ('material',mater, 'fes',bfes,...
            'integration_rule',eltyd(eix).surface_integration_rule));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = 0*geom; % zero out
        % Apply EBC's
        ebc_fenids=fenode_select (fens,struct('box',[0,0,-100*sR,100*sR,-100*sR,100*sR],'inflate',R/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+1;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[-100*sR,100*sR,0,0,-100*sR,100*sR],'inflate',R/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[-100*sR,100*sR,-100*sR,100*sR,0,0],'inflate',R/1000));
        ebc_prescribed=ones(1,length (ebc_fenids));
        ebc_comp=ebc_fenids*0+3;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        
        % Number equations
        u   = numberdofs (u);
        % Now comes the nonlinear solution
        tup = 1;
        u = u*0; % zero out the displacement
        utol = 1e-13*u.nfreedofs;
        
        gv=graphic_viewer;
        
        t=0; % time=load magnitude
        incr=1;
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
            femm1 =femm;
            while 1
                %                 prop = property_linel_iso (struct('E',E,'nu',nus(iter)));
                %                 mater = mater_defor_neohook(struct('property',prop));
                %                 femm1=set(femm,'mater',mater);
                %                 femm1=femm;                                           % Make a copy of the state
                Load=zeros(3,1); Load(3) =-Fmag*t;
                fi=force_intensity(struct('magn',Load));
                FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
                F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
                K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
                % Displacement increment
                du = scatter_sysvec(du, K\F);
                R0 = dot(F,gather_sysvec(du));
                F = FL + restoring_force(femm1,sysvec_assembler, geom,u1+du,u);       % Internal forces
                R1 = dot(F,gather_sysvec(du));
                a = R0/R1;
                if ( a<0 )
                    eta = a/2 +sqrt((a/2)^2 -a);
                else
                    eta =a/2;
                end
                eta=min( [eta, 1.0] );
                u1 = u1 + eta*du;   % increment displacement
                %                 draw(sfemm,gv, struct ('x', geom,'u',u1, 'facecolor','none'));
                %                 pause(1)
                disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
                if (norm(du) < utol) break; end;                    % convergence check
                iter=iter+1;
            end
            [ignore,femm] = restoring_force(femm,sysvec_assembler,geom,u1,u);        % final update
            disp(['    Converged for t=' num2str(t)]); % pause
            u = u1;                                               % update the displacement
            gv=reset(clear(gv,[]),[]);
            %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
            %             draw(sfemm,gv, struct ('x', geom,'u',0*u, 'facecolor','none'));
            fld = field_from_integration_points(femm, geom, u1, [], 'Cauchy',3);
            nvals=fld.values;%min(nvals),max(nvals)
            nvalsrange=[min(nvals),max(nvals)];
            dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
            draw(sfemm,gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink',1.0));
            draw(sfemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
            % draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
            colormap(cmap);
            cbh=colorbar;
            set(cbh,...
                'Position',[0.815 0.15 0.05 0.7],...
                'YLim',[0,1],...
                'YTick',[0,1],...
                'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
            set(get(cbh,'XLabel'),'String','\sigma_z');
            pause(0.5);
            incr = incr + 1;
        end
        clear K
        save ([mfilename '-nR' num2str(nR) '-nt' num2str(nt) '-' eltyd(eix).description])
    end

for eix = 1:length(eltyd)
    for i=1:length(nR)
        for j=1:length(nt)
            Compute(nR(i),nt(j));
        end
    end
end
end