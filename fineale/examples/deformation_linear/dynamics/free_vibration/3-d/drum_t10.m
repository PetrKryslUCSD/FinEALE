% Circular clamped plate, T10 elements.
E=0.1e6;% Pa
nu=0.3;
rho=1000;% kg
R= 25.0e-3;% m
t= 2.0e-3;% m
rand('state',0);% try to comment out this line and compare 
%                   results for several subsequent runs

% Mesh
[fens,fes] = T4_cylinderdel(t,R,1,5);
% Convert to quadratic tetrahedra
[fens,fes] = T4_to_T10(fens,fes);
% Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_linear_triax (struct('property',prop ));
    %     Make the finite element model machine
femm = femm_deformation_linear (struct ('material',mater, 'fes',fes,...
    'integration_rule',tet_rule (struct('npts',4))));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
for i=1:count(fens)
    if abs(R-norm(fens.xyz(i,2:3))) <0.005*R
        u   = set_ebc(u, i, [1], [], 0.0);
    end
end
u   = apply_ebc (u);
% Number equations
u   = numberdofs (u);
% Assemble the system matrix
K = stiffness(femm, sysmat_assembler_sparse, geom, u);
M = mass(femm, sysmat_assembler_sparse, geom, u);
%
neigvs = 4;
[W,Omega]=eigs(K,M,neigvs,'SM');
[Omegas,ix]=sort(diag(Omega));
Omega= diag(Omegas);

% 
lambda_ij2 = 10.22;
disp(['  Analytical eigenvalue 1: '  ...
    num2str(lambda_ij2/(2*pi*R^2)*sqrt (E*t^3/(12*rho*t*(1-nu^2)))) ' Hz']);
lambda_ij2 = 21.26;
disp(['  Analytical eigenvalue 2, 3: '  ...
    num2str(lambda_ij2/(2*pi*R^2)*sqrt (E*t^3/(12*rho*t*(1-nu^2)))) ' Hz']);
lambda_ij2 = 34.88;
disp(['  Analytical eigenvalue 4: '  ...
    num2str(lambda_ij2/(2*pi*R^2)*sqrt (E*t^3/(12*rho*t*(1-nu^2)))) ' Hz']);

% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
scale=0.065;
w=clone(u,'w'); % make a copy of u
bfes=mesh_boundary(fes);
for i=1:neigvs
    disp(['  Eigenvector ' num2str(i) ' frequency ' num2str(sqrt(Omega(i,i))/2/pi) ]);
    w = scatter_sysvec(w, W(:,ix(i)));
    wmag = magnitude(w);
    dcm=data_colormap(struct ('range',[min(wmag.values),max(wmag.values)], 'colormap',jet));
    colors=map_data(dcm, wmag.values);
    colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
    gv=reset (gv,[]);
    set(gca,'FontSize',16)
    camset(gv,[-0.1591    0.0085   -0.3252    0.0001    0.0043    0.0013   -0.8985    0.0237    0.4383    5.4322]);
    draw(bfes,gv, struct ('x', geom,'u',+scale*w,'colorfield',colorfield));
    text(0,1.1*R,1.1*R,['\omega_' num2str(i) '=' num2str(sqrt(Omega(i,i))/2/pi) ],'FontSize',24);
    axis off
%     saveas(gcf, ['drum_t4-' num2str(i) '.png'], 'png');
    pause(2); % next eigenvector
end

