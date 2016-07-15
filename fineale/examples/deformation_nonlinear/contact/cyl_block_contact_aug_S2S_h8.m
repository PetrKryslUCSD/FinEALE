% Rigid cylinder pushed into a block.
% Augmented Lagrangean contact approach.
function  pin_and_lug_contact_aug_S2S_h8
disp('Rigid cylinder pushed into a block. Augmented Lagrangean contact algorithm.');

c =clock; disp(['Started ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6)) ' '])
% Parameters:
U=physical_units_struct;
E=2*U.MEGA*U.PA;
nu=0.49;
L= 0.03*U.M; % in-plane dimension
W = 0.015*U.M; % in-plane dimension
a= 0.005*U.M; % cylinder radius
H = 0.01*U.M; % thickness of the block
nL=2*13;nH=10;nW=13; 
tolerance=min([L,W,a,H])/1e5;
tol = a*10e-7;
penalty = 1000.0*E;
utol=1e-9*U.M;
augtol=1e-12;
Rrat= 0.5;
surface_data.center=[0,0,-Rrat*a];
surface_data.direction=[0,0,+1];
surface_data.R2=a^2;

function  [penetration,normal] = get_penetration (surface_data,X,U)
    p=X-surface_data.center;
    P=p-surface_data.direction*dot(p,surface_data.direction);
    if (sum(P.*P)<=surface_data.R2)
        normal=(-1)*surface_data.direction;
        penetration =  - (p+U)*normal';
    else
        normal=[];% No penetration
        penetration =0;
    end
end

% Mesh
[fens,fes] = H8_block(L,W,H,nL,nW,nH);
fens.xyz=[fens.xyz(:,1)-L/2,fens.xyz(:,2),fens.xyz(:,3)-H];
%    gv=drawmesh({fens,fes},'fes', 'facecolor','r')

bdry_fes = mesh_boundary(fes,  []);
fixedfacelist =[fe_select(fens,bdry_fes,struct ('facing',true,'direction',[+1,0,0],'tolerance',tolerance)),...
    fe_select(fens,bdry_fes,struct ('facing',true,'direction',[-1,0,0],'tolerance',tolerance)),...
    fe_select(fens,bdry_fes,struct ('facing',true,'direction',[0,+1,0],'tolerance',tolerance)),...
    fe_select(fens,bdry_fes,struct ('facing',true,'direction',[0,0,-1],'tolerance',tolerance)),...
    ];
freefacelist =[fe_select(fens,bdry_fes,struct ('facing',true,'direction',[0,0,+1],'tolerance',tolerance)),...
    ];
symmfacelist =[fe_select(fens,bdry_fes,struct ('facing',true,'direction',[0,-1,0],'tolerance',tolerance)),...
    ];
penetration_fes=subset(bdry_fes,freefacelist);
        gv=drawmesh({fens,penetration_fes},'fes', 'facecolor','r')
    gv=drawmesh({fens,subset(bdry_fes,fixedfacelist)},'gv',gv,'fes', 'facecolor','y')
   gv=drawmesh({fens,subset(bdry_fes,symmfacelist)},'gv',gv,'fes', 'facecolor','b')




% Material
prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));
% Finite element block
femm = femm_deformation_linear_h8msgso(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
penetrationfemm = femm_deformation_linear_penetration_aug_lag (struct('fes',penetration_fes,...
    'integration_rule', trapezoidal_rule(struct('dim',2)),'get_penetration',@get_penetration,'surface_data',surface_data,'penalty',penalty));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = clone(geom,'u');
u   = u*0; % zero out
% Apply EBC's
ebc_fenids=connected_nodes(subset(bdry_fes,fixedfacelist));
u   = set_ebc(u, ebc_fenids, true, 1, 0.0);
u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
u   = set_ebc(u, ebc_fenids, true, 3, 0.0);
ebc_fenids=connected_nodes(subset(bdry_fes,symmfacelist));
u   = set_ebc(u, ebc_fenids, true, 2, 0.0);
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);


u   = u*0; % zero out

femm = associate_geometry(femm, geom);
tic;
Ks = stiffness(femm, sysmat_assembler_sparse, geom, u);
disp(['Stiffness matrix assembly ' num2str(toc)])
%     R = resultant_force(u,Fl)
du=0*u;
u1=0*u;
gv=reset (graphic_viewer,[]);
for augit= 1:25
    disp([' Augmented Iteration ' num2str(augit)])
    iter=1;
    while 1
        Fc = contact_loads(penetrationfemm, sysvec_assembler, geom, u1);
        %                         resultant_force(u,Fc)
        U1= gather_sysvec(u1);
        Fr  = -(Ks*U1);
        Kg =stiffness(penetrationfemm, sysmat_assembler_sparse, geom, u1);
        K = Ks +  Kg;
        du = scatter_sysvec(du, K\(Fc+Fr));
        u1 = u1 + du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(max(abs(du.values))) < utol) break; end;                    % convergence check
        iter=iter+1;
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
draw(bdry_fes,gv, struct ('x', geom, 'u', scale*u,'colorfield',colorfield));
draw_cylinder(gv, surface_data.center, surface_data.center+surface_data.direction*H, a, a,...
    struct('facecolor','b','tessel',40));
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on
end
