% Infinite plane with a circular hole.
%
% Infinite plane with a circular hole of radius a, subject to  
% uniform tension in the X direction. Analytical solution is available.
%
% This version spells out all the operations explicitly.
function plane_w_hole
    % Parameters:
    U=physical_units_struct;
    E=210000*U.MEGA*U.PA;
    nu=0.3;
    L= 0.3*U.M; % in-plane dimension
    W = 0.3*U.M; % in-plane dimension
    a= 0.15*U.M; % hole radius
    H = 0.01*U.M; % thickness of the plate
    nL=5;nH=1;nW=5;na=25;
    tol = a*10e-7;
    sigma0=1*U.MEGA*U.PA;
    graphics = true;

    % Mesh
    [fens,fes]=Q4_elliphole(a,a,L,W,nL,nW,na,[]);
    [fens,fes] = H8_extrude_Q4(fens,fes,nH,@(x,i)([x,0]+[0,0,H*i]));
    [fens,fes] = H8_to_H20(fens,fes);
    % Material
    prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
        mater = material_deformation_linear_triax (struct('property',prop));
        % Finite element block
    femm = femm_deformation_linear (struct ('material',mater, 'fes',fes,...
        'integration_rule', gauss_rule (struct('dim', 3,  'order', 2))));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = 0*geom; % zero out
    % Apply EBC's
    ebc_fenids=fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',tol));
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=ones(1,length (ebc_fenids))*1;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,0,0,-Inf,Inf],'inflate',tol));
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=ones(1,length (ebc_fenids))*2;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    ebc_fenids=[fenode_select(fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0],'inflate',tol)),fenode_select(fens,struct('box',[-Inf,Inf,-Inf,Inf,H,H],'inflate',tol))];
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=ones(1,length (ebc_fenids))*3;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    % Assemble the system matrix
    % Load
    bdry_fes = mesh_boundary(fes, []);
    bclx = fe_select(fens, bdry_fes, ...
        struct ('box',[L,L,-Inf,Inf,-Inf,Inf],'inflate',tol));
    lfemmx = femm_deformation_linear(struct ('mater',mater, 'fes',subset(bdry_fes,bclx),...
        'integration_rule',gauss_rule (struct('dim', 2,  'order', 2))));
    bcly = fe_select(fens, bdry_fes, ...
        struct ('box',[-Inf,Inf,W,W,-Inf,Inf],'inflate',tol));
    lfemmy = femm_deformation_linear(struct ('mater',mater, 'fes',subset(bdry_fes,bcly),...
        'integration_rule',gauss_rule (struct('dim', 2,  'order', 2))));
    % Assemble the system matrix
K =  stiffness(femm, sysmat_assembler_sparse, geom, u);
% Load
function sxx =sigmaxx(x)
r=norm(x(1:2));
th =atan2(x(2),x(1));
sxx = sigma0*(1-a^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*a^4/r^4*cos(4*th));
end
function syy =sigmayy(x)
r=norm(x(1:2));
th =atan2(x(2),x(1));
syy = -sigma0*(a^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*a^4/r^4*cos(4*th));
end
function sxy =sigmaxy(x)
r=norm(x(1:2));
th =atan2(x(2),x(1));
sxy = -sigma0*(a^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*a^4/r^4*sin(4*th));
end

fi=force_intensity(struct('magn',@(x) ([sigmaxx(x),sigmaxy(x),0])));
F = distrib_loads(lfemmx, sysvec_assembler, geom, u, fi, 2);
fi=force_intensity(struct('magn',@(x) ([sigmaxy(x),sigmayy(x),0])));
F = F + distrib_loads(lfemmy, sysvec_assembler, geom, u, fi, 2);
% Solve
u = scatter_sysvec(u, K\F);
% get(u,'values')

% Plot
gv=graphic_viewer;
gv=reset (gv,struct ([]));
scale=10000;
cmap = jet;
fld = field_from_integration_points(femm, geom, u, [], 'Cauchy',1);
nvals=fld.values/(U.MEGA*U.PA);%min(nvals),max(nvals)
nvalsrange=[min(nvals),max(nvals)];
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
draw(femm,gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink',1.0));
draw(femm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',0.8));
% draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red'));
colormap(cmap);
cbh=colorbar;
set(cbh,...
    'Position',[0.815 0.15 0.05 0.7],...
    'YLim',[0,1],...
    'YTick',[0,1],...
    'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
set(get(cbh,'XLabel'),'String','\sigma_x');
view (2)
set_graphics_defaults(gcf)

end