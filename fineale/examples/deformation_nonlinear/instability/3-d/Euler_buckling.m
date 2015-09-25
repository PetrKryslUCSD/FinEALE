function Euler_buckling
disp('Euler buckling: small strain nonlinear analysis');
% Parameters:
E=100;
nu=0.3;
L= 28; % Length of the beam
W = 1; % width
H = 3; % height
magn=0.026226627341330;
scale=5;
magn_cr = E*min([H,W])^3*max([H,W])/12*pi^2/(2*L)^2/W/H

neigvs = 4;

% Mesh
% [fens,gcells] = block(W, H, L, 2,2, 18);
% [fens,gcells] = block(W, H, L, 2,2, 8);
[fens,gcells] = H8_block(W, H, L, 1,1, 12);
% [fens,gcells] = block(W, H, L, 1,1,2);
[fens,gcells] = H8_to_H27(fens,gcells);
% Material
prop = property_linel_iso(struct('E',E,'nu',nu));
mater = mater_defor_ss_linel_triax(struct('property',prop));
% Finite element block
feb = feblock_defor_ss_nonlinear (struct ('mater',mater, 'gcells',gcells,...
    'integration_rule',gauss_rule(3,3)));
% Geometry
geom = field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out

% Apply EBC's
ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0,],'inflate',1/1000));
ebc_prescribed=ones(1,length (ebc_fenids));
ebc_comp=[];
ebc_val=ebc_fenids*0;
u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
u   = apply_ebc (u);
% Number equations
u   = numbereqns (u);
% Assemble the system matrix
ems = stiffness(feb, geom, u, u);
K = start (sparse_sysmat, get(u, 'neqns'));
K = assemble (K, ems);
% Load
bdry_gcells = mesh_bdry(gcells, []);
bcl = gcell_select(fens, bdry_gcells, ...
    struct ('box',[-Inf,Inf,-Inf,Inf,L L],'inflate',H/1000));
lfeb = feblock_defor_ss(struct ('mater',mater, 'gcells',subset(bdry_gcells,bcl),...
    'integration_rule',gauss_rule(2,2)));
fi= force_intensity(struct ('magn',[0;0; -magn]));
F = start (sysvec, get(u, 'neqns'));
F = assemble (F, distrib_loads(lfeb, geom, u, fi,2));
% Solve
u = scatter_sysvec(u, K\F);


ulin = u; %  linear displacement solution

disp('Stability analysis');
[feb1,evs] = restoring_force(feb,geom,ulin,0*ulin);       % Internal forces
ems = stiffness_geo(feb1, geom, ulin, ulin);            % Geometric stiffness
K_geo = start (sparse_sysmat, get(u, 'neqns'));
K_geo = assemble (K_geo, ems);

options.tol =1e-26;
options.maxit = 5000;
options.disp = 0;
[W,Omega]=eigs(-(get(K_geo,'mat')+get(K_geo,'mat')')/2,(get(K,'mat')+get(K,'mat')')/2,...
    neigvs,'lm', options);
diag(Omega);
%         [Omegas,ix]=sort(1./diag(Omega)) ;
Omegas=(1./(diag(Omega)));


gv=graphic_viewer;
cmap = jet(16);
w=clone(ulin,'w');
for i=1:neigvs
    disp(['  Eigenvector ' num2str(i) ' eigenvalue ' num2str(Omegas(i)) ]);
    gv=reset(clear(gv,[]),[]);
            
    w = scatter_sysvec(w, W(:,i));
    v = get(w,'values');
    wmag = field(struct ('name',['wmag'], 'dim', 1, 'nfens',count(fens)));
    wmag = scatter(wmag,1:get(wmag,'nfens'),sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2));
    nvals =gather(wmag,1:get(wmag,'nfens'),'values');
    nvalsrange=[min(nvals),max(nvals)];
    dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
    colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
    draw(feb,gv, struct ('x', geom, 'u', scale*w,'colorfield',colorfield, 'shrink',1.0));

    pause(1); % next eigenvector
end
o=somel(sort(Omegas),1:4);
oref= [ 1.009374713026403
   8.981963216219038
   9.059348763440447
  25.027908601997087];
    assignin('caller','faesor_test_passed',(norm(o-oref)<1e-3));
end