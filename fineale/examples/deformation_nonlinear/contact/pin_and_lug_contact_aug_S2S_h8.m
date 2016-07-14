% Lug and pin.
% Augmented Lagrangean contact approach.
function  pin_and_lug_contact_aug_S2S_h8
disp('Lug and pin: the cylindrical pin surface is  rigid. Augmented Lagrangean contact algorithm.');

c =clock; disp(['Started ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
% Parameters:
U=physical_units_struct;
E=210000*U.MEGA*U.PA;
nu=0.3;
L= 0.03*U.M; % in-plane dimension
W = 0.03*U.M; % in-plane dimension
a= 0.015*U.M; % hole radius
H = 0.005*U.M; % thickness of the plate
nL=5;nH=4;nW=5;na=15;tolerance=min([L,W,a,H])/1e5;
tol = a*10e-7;
sigma0=100*U.MEGA*U.PA;
sscale=1.0;
penalty = 0.1*E;
utol=1e-9*U.M;
augtol=1e-12;
Rrat= 0.9599;
surface_data.center=[(Rrat-1.0)*a,0,0];
surface_data.R=Rrat*a;

    function  [penetration,normal] = get_penetration (surface_data,X,U)
        normal = get_normal(surface_data,X);
        p=X-surface_data.center;
        penetration = surface_data.R - dot(p,normal) - dot(U,normal);
    end
    function normal = get_normal(surface_data,x)
        x(3)=0;
        normal=(x-surface_data.center);
        normal= normal/norm(normal);
    end

% Mesh
[fens,fes]=Q4_elliphole(a,a,L,W,nL,nW,na,[]);
[fens,fes] = H8_extrude_Q4(fens,fes,nH,@(x,i)([x,0]+[0,0,H/nH*i]));
[fens1,fes1] = mirror_mesh(fens, fes,...
    [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
[fens,fes1,fes2] = merge_meshes(fens, fes, fens1,    fes1, tolerance);
fes=cat(fes1,fes2);

bdry_fes = mesh_boundary(fes,  []);
pullfacelist =fe_select(fens,bdry_fes,struct ('box',[L L -inf inf -inf inf],'inflate',tolerance));
holefacelist =fe_select(fens,bdry_fes,struct ('cylinder',[0 0 0 a 20*H],'inflate',tolerance));
clampfacelist =fe_select(fens,bdry_fes,struct ('box',[-L -L -inf inf -inf inf],'inflate',tolerance));
penetration_fes=subset(bdry_fes,holefacelist);
cn = connected_nodes(subset(bdry_fes,clampfacelist));
Spring_fes=fe_set_P1 (struct('conn',cn));
cfnli=fenode_select(fens,struct ('box',[-a -a 0 0 0 0],'inflate',tolerance));
%         gv=drawmesh({fens,bdry_fes},'fes', 'facecolor','r')
%     gv=drawmesh({fens,subset(bdry_fes,holefacelist)},'gv',gv,'fes', 'facecolor','y')
%  draw_cylinder(gv, surface_data.center-[0,0,2*H], surface_data.center+[0,0,2*H], surface_data.R, surface_data.R,...
%          struct('facecolor','b','tessel',40));



% Material
prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));
% Finite element block
femm = femm_deformation_linear_h8msgso(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
penetrationfemm = femm_deformation_linear_penetration_aug_lag (struct('fes',penetration_fes,...
    'integration_rule', trapezoidal_rule(struct('dim',2)),'get_penetration',@get_penetration,'surface_data',surface_data,'penalty',penalty));
Springfemm = femm_deformation_spring_grounded (struct('fes',Spring_fes,...
    'translation_stiffness_matrix',penalty/1e9*diag([1,0,0])));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[-inf inf -inf inf 0 0],'inflate',tolerance));
u   = set_ebc(u, ebc_fenids, true, 3, 0.0);
ebc_fenids=fenode_select (fens,struct('box',[-inf inf 0 0 -inf inf],'inflate',tolerance));
u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
%     ebc_fenids=fenode_select (fens,struct('box',[-a -a 0 0 0 0],'inflate',tolerance));
%     u   = set_ebc(u, ebc_fenids, true, 1, 0.0);
%     ebc_fenids=connected_nodes(subset(bdry_fes,clampfacelist));
%     u   = set_ebc(u, ebc_fenids, true, [], 0.0);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);


u   = u*0; % zero out

femm = associate_geometry(femm, geom);
tic;
Ks = stiffness(femm, sysmat_assembler_sparse, geom, u);
Ksp =stiffness(Springfemm, sysmat_assembler_sparse, geom, u);
disp(['Stiffness matrix assembly ' num2str(toc)])
% Load
efemm = femm_deformation_linear (struct ('material',mater, ...
    'fes',subset(bdry_fes,pullfacelist),...
    'integration_rule',gauss_rule (struct('dim',2,'order', 2))));
fi=force_intensity(struct('magn',[sigma0;0;0]));
Fl = distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
%     R = resultant_force(u,Fl)
du=0*u;
u1=0*u;
gv=reset (graphic_viewer,[]);
for augit= 1:5
    disp([' Augmented Iteration ' num2str(augit)])
    iter=1;
    while 1
        Fc = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
        %                         resultant_force(u,Fc)
        U1= gather_sysvec(u1);
        Fr  = -(Ksp*U1+Ks*U1);
        Kg =stiffness(penetrationfemm, sysmat_assembler_sparse, geom, u1);
        K = Ks + Ksp + Kg;
        du = scatter_sysvec(du, K\(Fc+Fr+Fl));
        u1 = u1 + du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(max(abs(du.values))) < utol) break; end;                    % convergence check
        iter=iter+1;
        u1.values(cfnli,:)
    end
    prevlm=penetrationfemm.lm;
    [~,penetrationfemm] = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
    currlm=penetrationfemm.lm;
    disp(['   Augmented Lagrange It. ' num2str(augit) ': ||dlm||/||lm||=' num2str(norm(currlm-prevlm)/norm(currlm))]);
    if (norm(currlm-prevlm)<=augtol*norm(currlm)), break,end
    
    
    %     gv=reset (graphic_viewer,struct('limits',inflate_box(bounding_box(fens.xyz),a)));
    %     draw(bdry_fes,gv, struct ('x', geom, 'u', sscale*u1,'facecolor','none'));
    %     draw(subset(bdry_fes,holefacelist),gv, struct ('x', geom, 'u', (0)*u1,'facecolor','y'));
    %     LMForces=0*u1; LMForces=(1/700)*scatter_sysvec(LMForces,lm);
    %     show_field_as_arrows(gv,struct('x',geom, 'u', LMForces))
    %     ApplForces=0*u1; ApplForces=(1/700)*scatter_sysvec(ApplForces,Fl);
    %     show_field_as_arrows(gv,struct('x',geom, 'u', ApplForces))
    %     pause(1);
end

u=u1;% Now we have our final solution
%         if isempty(progressbar((step/nsteps),h)), return, end;

c =clock; disp(['Stopped ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
scale=1;
%         draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','none'));
w = u;
w = w*(1/norm(w ));
wv =w.values;
wmag = real (sqrt(wv(:,1).^2+wv(:,2).^2+wv(:,3).^2));
dcm=data_colormap(struct ('range',[min(wmag),max(wmag)], 'colormap',jet));
colors=map_data(dcm, wmag);
colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
draw(bdry_fes,gv, struct ('x', geom, 'u', sscale*u,'colorfield',colorfield));
draw_cylinder(gv, surface_data.center-[0,0,2*H], surface_data.center+[0,0,2*H], surface_data.R, surface_data.R,...
    struct('facecolor','b','tessel',40));
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on
end
