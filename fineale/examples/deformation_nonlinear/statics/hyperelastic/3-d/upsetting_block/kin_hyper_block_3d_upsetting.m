function kin_hyper_block_3d_upsetting(ntmult)
if ~exist('ntmult','var'), ntmult=5; end
% A new locking-free brick element technique for large
% deformation problems in elasticity p
% S. Reese a, *, P. Wriggers a , B.D. Reddy
%     kappa=400889.8;
E = 210000;
nu=0.4999;;

umag= 0.5;
gtol=1e-6;
utol = 1e-12;
maxiter=12;
tup = 1; dt= 0.01;
graphics = ~false;
scale=1;
nt = [ntmult*[2,1,2]];
%     nt = [4*[2,2,2]];

prop = property_deformation_neohookean(struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));

    function [fens,fes] = h8H8(nt)
        [fens,fes]=H8_block(2.0, 1.0, 1.0, nt(1),nt(2),nt(3));
    end


eix=1;
clear eltyd

eltyd(eix).description='H8MSGSO(U)';
eltyd(eix).mf =@h8H8;
eltyd(eix).blf =@(fes)femm_deformation_nonlinear_h8msgso(struct ('material',mater, ...
    'fes',fes,  'match_stabilization',true,...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
eix=eix+1;


    function [neqns,energy] =Compute(nt)
        %          Mesh
        [fens,fes] = eltyd(eix).mf (nt);
        %         [fens, fes]=renum_mesh(fens, fesAnd, 'symrcm');
        femm = eltyd(eix).blf (fes);
        bfes =mesh_boundary(fes, []);
        sfemm = femm_deformation (struct ('material',mater, 'fes',bfes,...
            'integration_rule',eltyd(eix).surface_integration_rule));
        %         [fens, fes]=renum_mesh_renumber(fens, fes,...
        %             renum_mesh_numbering(fens,genconn(feb), 'symrcm'));
        %         [fens, fes]=renum_mesh_renumber(renum_mesh_numbering(fens,genconn(feb), 'symamd'),fens, fes);
        %         icl = fe_select(fens, bfes, struct('box', [0,1,0,1, 1,  1],'inflate',1/100));
        cncl = fenode_select(fens, struct('box', [0,0,0,0,1,1],'inflate',gtol));
        %          gv=drawmesh( {fens,fes},'facecolor','none'); gv=drawmesh( {fens,bfes(icl)},'gv',gv,'facecolor','blue');
        %         view(3); pause(1)
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = 0*geom; % zero out
        % Apply EBC's
        %         Bottom surface: Clamped
        ebc_fenids=fenode_select (fens,struct('box',[-inf,inf,-inf,inf,0,0],'inflate',gtol));
        u   = set_ebc(u, ebc_fenids, true, [], 0.0);
        %         Side surface
        ebc_fenids=fenode_select (fens,struct('box',[0,0,-inf,inf,-inf,inf],'inflate',gtol));
        u   = set_ebc(u, ebc_fenids, true,  1, 0.0);
        ebc_fenids=fenode_select (fens,struct('box',[2,2,-inf,inf,-inf,inf],'inflate',gtol));
        u   = set_ebc(u, ebc_fenids, true,  1, 0.0);
        % Top: moving punch surface
        top_ebc_fenids=fenode_select (fens,struct('box',[0,1,-inf,inf,1,1],'inflate',gtol));
        u   = set_ebc(u, top_ebc_fenids, true, 3, 0.0);
        %Apply EBC
        u   = apply_ebc (u);
        
        % Number equations
        u   = numberdofs (u);
        % Now comes the nonlinear solution
        u = u*0; % zero out the displacement
        utol =         utol*u.nfreedofs;
        
        if (graphics),
            gv=reset(clear(graphic_viewer,[]),[]);
            cmap = jet;
            peekb  = uicontrol('Parent',gcf,'Style','pushbutton',...
                'Units','Normalized','Position',[.9 .0 .1 .05],...
                'String','Peek','Value',0,...
                'TooltipString','Invoke command line in order to peek at data',...
                'Callback',@(h,varargin)keyboard);
        end
        
        t=0; % time=load magnitude
        ts={t};
        us={u};
        femm  =associate_geometry(femm,geom);
        while (t <= tup)
            disp(['Time ' num2str(t) ]); % pause
            % Initialization
            u1 = u; % guess
            u1   = set_ebc(u1, top_ebc_fenids, true, 3, -umag*(t+dt)/tup);
            u1 = apply_ebc(u1);
            du = 0*u; % this will hold displacement increment
            du = apply_ebc(du);
            
            iter=1;
            converged=false;
            while ~converged
                FL=zeros(u.nfreedofs,1);
                F = FL + restoring_force(femm,sysvec_assembler, geom,u1,u);       % Internal forces
                K = stiffness(femm, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm, sysmat_assembler_sparse, geom, u1,u);
                % Displacement increment
                du = scatter_sysvec(du, K\F);
                eta= 1.0;
                R0 = dot(F,gather_sysvec(du));
                F = FL + restoring_force(femm,sysvec_assembler, geom,u1+du,u);       % Internal forces
                R1 = dot(F,gather_sysvec(du));
                a = R0/R1;
                if ( a<0 )
                    eta = a/2 +sqrt((a/2)^2 -a);
                else
                    eta =a/2;
                end
                eta=real(eta);
                if (imag(eta)~=0)
                    disp('######################  Inverted elements?')
                end
                eta=min( [eta, 1.0] );
                u1 = u1 + eta*du;   % increment displacement
                if graphics, pause(0.5); Cam =camget(gv); end
                disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
                %                 draw(sfemm,gv, struct ('x', geom, 'u', u1,'facecolor','none'));
                if (max(abs(du.values)) < utol)
                    [~,femm]  =restoring_force(femm,sysvec_assembler,geom,u1,u);
                    disp(['    Converged for t=' num2str(t)]); % pause
                    u = u1;       % update the displacement
                    dt=dt*maxiter/iter;
                    if graphics
                        gv=reset(clear(gv,[]),[]);
                        camset (gv,Cam);
                        if (0)
                            fld = field_from_integration_points(femm, geom, u1, [], 'pressure',1);
                            nvals=fld.values;%min(nvals),max(nvals)
                            nvalsrange=[min(nvals),max(nvals)];
                            dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
                            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
                            draw(sfemm,gv, struct ('x', geom, 'u', u1,'colorfield',colorfield, 'edgecolor',0.6*[1,1,1],'linewidth',2));
                        else
                            draw(sfemm,gv, struct ('x', geom, 'u', u,'facecolor','red'));
                        end
                        pause(0.5); Cam =camget(gv);
                    end
                    break;
                end;                    % convergence check
                iter=iter+1;
                if (iter >maxiter)% Restart the iteration with a shorter step
                    disp(['    Restart for t=', num2str(t), ' dt=',num2str(dt)]); % pause
                    dt=dt/2;
                    u1 =u; % restart
                    u1   = set_ebc(u, top_ebc_fenids, true, 3, -umag*(t+dt)/tup);
                    u1 = apply_ebc(u1);
                    iter=1;
                end
                
            end
            
            t=t+dt;
            ts{end+1} =t;
            us{end+1} =u;
            utip=gather_values(u,cncl);
            disp(['    Center compression =' num2str(-utip(3))]); % pause
        end % Increment loop
        
                
        clear K
        save ([mfilename  '-nt' num2str(nt) '-' eltyd(eix).description])
    end

for eix = 1:length(eltyd)
    Compute(nt);
end
end