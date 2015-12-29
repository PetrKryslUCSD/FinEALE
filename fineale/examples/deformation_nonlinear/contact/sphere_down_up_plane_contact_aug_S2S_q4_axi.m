% Fear squashed between two rigid planes
% Penalty contact approach.
function  Sphere_plane_contact_aug_S2S_q4_axi
disp('Sphere squashed between two rigid planes. Augmented Lagrangean contact algorithm.');

c =clock; disp(['Started ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
% Parameters:
U=physical_units_struct;
E=210000*U.MEGA*U.PA;
nu=0.3;
a= 0.015*U.M; % sphere radius
na=75;tolerance=min([a])/1e5;
sscale=1.0;
penalty = 500*E;
utol=1e-15*U.M;
augtol=1e-12;
surface_data.center=[0,0.9*a];

    function  [penetration,normal] = get_penetration (surface_data,X,U)
        normal = [0.0,-1.0];
        p=X-surface_data.center;
        penetration = - dot(p,normal) - dot(U,normal);
    end

% Mesh
[fens,fes]=Q4_circle_n(a,na,[]);
fes=fe_set_Q4(struct('conn',fes.conn,'axisymm', true));
bdry_fes = mesh_boundary(fes, struct('axisymm', true));
xfacelist =fe_select(fens,bdry_fes,struct ('box',[0 0  -inf inf],'inflate',tolerance));
yfacelist =fe_select(fens,bdry_fes,struct ('box',[-inf inf  0 0],'inflate',tolerance));
sfacelist=setdiff(1:count(bdry_fes),[xfacelist,yfacelist ]);
penetration_fes=subset(bdry_fes,sfacelist);
%         gv=drawmesh({fens,fes},'fes', 'facecolor','r')
%         view(2); labels
%     gv=drawmesh({fens,subset(bdry_fes,holefacelist)},'gv',gv,'fes', 'facecolor','y')
%  draw_cylinder(gv, surface_data.center-[0,0,2*H], surface_data.center+[0,0,2*H], surface_data.R, surface_data.R,...
%          struct('facecolor','b','tessel',40));



% Material
prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_biax (struct('property',prop, 'reduction','axisymm' ));
% Finite element block
femm = femm_deformation_linear(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',2, 'order',2))));
penetrationfemm = femm_deformation_linear_penetration_aug_lag (struct('fes',penetration_fes,...
    'integration_rule', gauss_rule(struct('dim',1, 'order',4)),'get_penetration',@get_penetration,'surface_data',surface_data,'penalty',penalty));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=connected_nodes(subset(bdry_fes,xfacelist));
u   = set_ebc(u, ebc_fenids, true, 1, 0.0);
ebc_fenids=connected_nodes(subset(bdry_fes,yfacelist));
u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);


u   = u*0; % zero out

tic;
Ks = stiffness(femm, sysmat_assembler_sparse, geom, u);
disp(['Stiffness matrix assembly ' num2str(toc)])

du=0*u;
u1=0*u;
% gv=reset (graphic_viewer,[]);
for augit= 1:5
    disp([' Augmented Iteration ' num2str(augit)])
    iter=1;
    while 1
        Fc = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
        U1= gather_sysvec(u1);
        Fr  = -(Ks*U1);
        Kg =stiffness(penetrationfemm, sysmat_assembler_sparse, geom, u1);
        K = Ks +  Kg;
        du = scatter_sysvec(du, K\(Fc+Fr));
        u1 = u1 + du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(max(abs(du.values))) < utol) break; end;                    % convergence check
        iter=iter+1;
        
        %     gv=reset (graphic_viewer,struct('limits',inflate_box(bounding_box(fens.xyz),a)));
        %     draw(bdry_fes,gv, struct ('x', geom, 'u', u1,'facecolor','none'));
        %     draw(subset(bdry_fes,holefacelist),gv, struct ('x', geom, 'u', (0)*u1,'facecolor','y'));
        %     glForces=0*u1; glForces=(1/7000)*scatter_sysvec(glForces,gl);
        %     show_field_as_arrows(gv,struct('x',geom, 'u', glForces))
        %     ApplForces=0*u1; ApplForces=(1/7000)*scatter_sysvec(ApplForces,Fl);
        %     show_field_as_arrows(gv,struct('x',geom, 'u', ApplForces))
        %     labels; pause(1);
        
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

% Plot
gv=graphic_viewer;
gv=reset (gv,struct('limits', bounding_box(geom.values)));
scale=1;
%         draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','none'));
w = u;
wv =w.values;
wmag = real (sqrt(wv(:,1).^2+wv(:,2).^2));
dcm=data_colormap(struct ('range',[min(wmag),max(wmag)], 'colormap',parula));
colors=map_data(dcm, wmag);
colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
draw(fes,gv, struct ('x', geom, 'u', sscale*u,'colorfield',colorfield));
draw_polyline(gv, [surface_data.center;surface_data.center+[a,0]], [1,2], struct('linewidth',3))
draw_colorbar(gv,struct('range',dcm.range, 'colormap',dcm.colormap,'position',[0.751,0.2,0.05,0.7],'label','umag'))
xlabel('X');
ylabel('Y');
hold on
pause(1)

surface_data.center=[0,1.001*a];
penetrationfemm.surface_data=surface_data;

du=0*u;
% gv=reset (graphic_viewer,[]);
for augit= 1:5
    disp([' Augmented Iteration ' num2str(augit)])
    iter=1;
    while 1
        Fc = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
        U1= gather_sysvec(u1);
        Fr  = -(Ks*U1);
        Kg =stiffness(penetrationfemm, sysmat_assembler_sparse, geom, u1);
        K = Ks +  Kg;
        du = scatter_sysvec(du, K\(Fc+Fr));
        u1 = u1 + du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(max(abs(du.values))) < utol) break; end;                    % convergence check
        iter=iter+1;
        
        %     gv=reset (graphic_viewer,struct('limits',inflate_box(bounding_box(fens.xyz),a)));
        %     draw(bdry_fes,gv, struct ('x', geom, 'u', u1,'facecolor','none'));
        %     draw(subset(bdry_fes,holefacelist),gv, struct ('x', geom, 'u', (0)*u1,'facecolor','y'));
        %     glForces=0*u1; glForces=(1/7000)*scatter_sysvec(glForces,gl);
        %     show_field_as_arrows(gv,struct('x',geom, 'u', glForces))
        %     ApplForces=0*u1; ApplForces=(1/7000)*scatter_sysvec(ApplForces,Fl);
        %     show_field_as_arrows(gv,struct('x',geom, 'u', ApplForces))
        %     labels; pause(1);
        
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

% Plot
gv=graphic_viewer;
gv=reset (gv,struct('limits', bounding_box(geom.values)));
scale=1;
%         draw(femm,gv, struct ('x', geom, 'u', 0*u, 'facecolor','none'));
w = u;
wv =w.values;
wmag = real (sqrt(wv(:,1).^2+wv(:,2).^2));
dcm=data_colormap(struct ('range',[min(wmag),max(wmag)], 'colormap',parula));
colors=map_data(dcm, wmag);
colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
draw(fes,gv, struct ('x', geom, 'u', sscale*u,'colorfield',colorfield));
draw_polyline(gv, [surface_data.center;surface_data.center+[a,0]], [1,2], struct('linewidth',3))
draw_colorbar(gv,struct('range',dcm.range, 'colormap',dcm.colormap,'position',[0.751,0.2,0.05,0.7],'label','umag'))
xlabel('X');
ylabel('Y');
set_presentation_defaults
hold on
end
