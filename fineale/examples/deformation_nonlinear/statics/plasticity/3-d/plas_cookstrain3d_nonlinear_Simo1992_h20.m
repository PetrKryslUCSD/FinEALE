function plas_cookstrain3d_nonlinear_Simo1992
% Finite-deformation plastic  Cook trapezoidal panel. Saturation-hardening plasticity.
%
% Reference:
% Simo, J. C. and F. Armero (1992). "GEOMETRICALLY NONLINEAR ENHANCED
% STRAIN MIXED METHODS AND THE METHOD OF INCOMPATIBLE MODES." International
% Journal for Numerical Methods in Engineering 33(7): 1413-1449.  
%     2.0   3.8315
%     4.0    6.4359
%     8.0    6.7527
%    16.0539    6.9073
%    32.0    6.9227
% H8MSGSO
% 2.0  7.5731
% 4.0   6.0748
% 8.0 6.9232
% 16.0   6.8500
% 32.0   6.9334

graphics= ~true;
scale=1;
epscale=0.0002*scale;
doeig=~true;
aspect=1;
smult=0;
ne= [1,2];
nu= 0.29;
E=206.900;%MPa
sigma_y=0.45; sigma_res=0.715;
delta_expon=16.93; K_lin=0.12924; beta=1;
magn =5/16;
nincr =10;
maxdu_tol = 48/10000;
tol=48/10000;
%
convutip=6.94;


%  Create the mesh and initialize the geometry
n=32;
[fens,fes] = H20_block(48,44,16/n, n, aspect*n, 1);
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


% Package model data
clear model_data;
model_data.fens =fens;

clear region

prop = property_deformation_plasticity_saturation_hardening(struct('E',E,'nu',nu,'sigma_y',sigma_y,'sigma_res',sigma_res,'delta_expon',delta_expon,'K_lin',K_lin,'beta',beta));
mater = material_deformation_unrotated_j2_hard_triax(struct('property',prop));

region.femm= femm_deformation_nonlinear(struct ('material',mater,'fes',fes,'integration_rule',gauss_rule(struct('dim',3,'order',2))));
surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
model_data.region{1} =region;

%  Clamped cross-section
clear essential
essential.component= [1,2,3];
essential.fixed_value= 0;
essential.node_list = fenode_select (fens,struct ('box',[0,0,0,44,0,max(fens.xyz(:,3))],'inflate',44/10000));
model_data.boundary_conditions.essential{1} = essential;

clear essential
essential.component= [3];
essential.fixed_value= 0;
essential.node_list = (1:count(fens));
model_data.boundary_conditions.essential{2} = essential;

clear traction
bdry_fes = mesh_boundary(fes, []);
bcl = fe_select(fens, bdry_fes,struct ('box',[48,48,44,60,0,max(fens.xyz(:,3))],'inflate',48/1000));
traction.fes =subset(bdry_fes,bcl);;
traction.traction= [0;magn;0];
traction.integration_rule =surface_integration_rule;
model_data.boundary_conditions.traction{1} = traction;

% If online graphics  is needed, initialize some variables
if (graphics),
    bdry_fes = mesh_boundary(fes, []);
    sfemm = femm_deformation (struct ('material',[], 'fes',bdry_fes,...
        'integration_rule',[]));
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    %Cam= 1.0e+03 *[-0.9065   -1.3161    1.0356    0.1802    0.1000    0.0050         0         0    0.0010    0.0078];
end

%model_data.result = fenode_select (fens,struct('box',[0,0,R,R,0,2*W],'inflate',W/1000));
utip = fenode_select (fens,struct('box',[48,48,60,60,0,1],'inflate',2/100));
% Select the solver options
model_data.load_multipliers=(1:nincr)/nincr*1.0;
model_data.maxdu_tol  =maxdu_tol;;
model_data.line_search  = true;
model_data.iteration_observer =@iteration_observer;
us={}; Uy=[];
model_data.load_increment_observer =@load_increment_observer;
% Call the nonlinear deformation solver
model_data =deformation_nonlinear_statics(model_data);

%-----Info for VM plot------------------------------
%         femm=model_data.region{1}.femm;
%         geom=model_data.geom;
%         u1=model_data.u;
%         scale=1;
%         save ([mfilename  '-T10-8-2-18-Hardening'])
%----------------------------------------------------

%     Report results
%     Center_fenids=fenode_select (fens,struct('box',[L,L,W/2, W/2,-Inf,Inf],'inflate',1/1000));
%     u1s=[];
%     for j=1:length(us)
%         u1s=[u1s,mean(gather_values(us{j},enl))];
%     end
%     u1s  =reshape(u1s',[3,nincr])';

% Observer function to be called when convergence is reached.
    function load_increment_observer(lambda,model_data)
        fprintf(1,'lambda=%g\n',lambda);
        if graphics
            gv=reset(clear(gv,[]),[]);
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', 0*model_data.un1,'facecolor','none', 'shrink',1.0));
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', model_data.un1,'facecolor','y', 'shrink',1.0));
            %             camset (gv,Cam);
            interact(gv);
            pause(0.5); Cam =camget(gv);
        end
        Uy=[ Uy,mean(model_data.un1.values(utip,2))]
        
        
        us{end+1} =model_data.un1;
    end

% Iteration of observer can be called as the solution is being computed.
    function iteration_observer(lambda,iter,du,model_data)
        fprintf(1,'%d: %g\n',iter,norm(du));
        %         if 1 && graphics
        %             gv=reset(clear(gv,[]),[]);
        %             draw(sfemm,gv, struct ('x', model_data.geom, 'u', 0*model_data.u,'facecolor','none', 'shrink',1.0));
        %             draw(sfemm,gv, struct ('x', model_data.geom, 'u', scale*model_data.u,'facecolor','y', 'shrink',1.0));
        %             camset (gv,Cam);
        %             interact(gv);
        %             pause(0.5); Cam =camget(gv);
        %         end
        if (graphics)
            id.comp= 1;
            id.container=-Inf;
            id=inspect_integration_points(model_data.region{1}.femm, model_data.geom, model_data.un1, model_data.un, model_data.dt, [],...
                (1:count (fes)), struct ('output',['equiv_pl_def']),...
                @mx,id);
            max_equiv_pl_def=id.container;
            id.container=Inf;
            id=inspect_integration_points(model_data.region{1}.femm, model_data.geom, model_data.un1, model_data.un, model_data.dt, [], ...
                (1:count (fes)), struct ('output',['equiv_pl_def']),...
                @mn,id);
            min_equiv_pl_def =id.container;
            dcm=data_colormap(struct ('range',[min_equiv_pl_def,max_equiv_pl_def], 'colormap',jet));
            gv=reset(clear(gv,[]),[]);
            title (['Iteration ' num2str(iter), ', max equiv pl def=',num2str(max_equiv_pl_def) ])
            %                 camset (gv,1.0e+002 *[ -2.1416   -1.4296    3.3375    0.1981    0.1191   -0.0063    0.0006    0.0004    0.0006 0.0039]);
            draw(model_data.region{1}.femm,gv, struct ('x', model_data.geom,...
                'u',scale*model_data.un1, 'facecolor','none'));
            draw_integration_points(model_data.region{1}.femm,gv,struct ('x',model_data.geom,...
                'un1',model_data.un1,'un',model_data.un,'dt',model_data.dt,'u_scale',scale, 'scale',epscale,'output',['equiv_pl_def'],'component',1,'data_cmap', dcm));
            drawnow;
            pause(0.1)
        end
        
        function id= mn(id,out,xyz,U,pc)
            id.container=min(out(id.comp), id.container);
        end
        
        function id= mx(id,out,xyz,~,pc)
            id.container=max(out(id.comp), id.container);
        end
    end
end
