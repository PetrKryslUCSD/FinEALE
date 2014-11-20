% Buckling of a clamped box at the bottom loaded uniformly other top.
% Reference solution from ABAQUS 38.366,30.821,30.824,23.489
% with about 45,000 degrees of freedom
function Box_buckling_8
    disp('Linear box buckling analysis');
    % Parameters:
    E=210e9;
    nu=0.3;
    Thickness= 0.1; iDim = 10/2-Thickness/2; eDim = iDim +Thickness; Length= iDim+eDim;
    
    nD= 2;nt=1;nL=4*nD;
    %   nD= 4;nt=1;nL=2*nD;
%     nD= 6;nt=1;nL=2*nD;
%     nD= 8;nt=1;nL=2*nD;
    
    D=E*Thickness^3/12/(1-nu^2);
    magn=D/Thickness/(iDim+eDim)^2;
    scale=5000;
    neigvs = 4;
    uscale=20000;
    Reference= [23.489,30.821,30.824,38.366,];
    
    % Mesh
    
    % [fens,fes]=block(rA, Thickness, Length, nc,nt,na);
    [fens,fes] = Q4_quadrilateral([iDim,0;eDim,0;eDim,eDim;iDim,iDim],nt,nD,[]);
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,1], [0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, Thickness/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [0,1], [0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, Thickness/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = mirror_mesh(fens, fes, [1,0], [0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, Thickness/1000);
    fes=cat(fes1,fes2);
    function xyz= up(xyz, layer)
        xyz= [xyz, +(layer/nL)*Length];
    end
[fens,fes] = Q4_refine(fens,fes);
    [fens,fes] = H8_extrude_Q4(fens,fes,nL,@up);
    
    drawmesh({fens,fes},'fes')
    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',1));
        mater = material_deformation_stvk_triax (struct('property',prop ));
        % Finite element block
    femm = femm_deformation_nonlinear_h8msgso (struct ('material',mater, 'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',3,'order',2))));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = 0*geom; % zero out
    
    % Apply EBC's
    ebc_fenids=fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0,],'inflate',Thickness/1000));
    ebc_prescribed=ones(1,length (ebc_fenids));
    ebc_comp=[];
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    disp([num2str(u.nfreedofs) ' degrees of freedom'])
    
    % Update the FEMM
    femm  =update(femm,geom,u,u);
    % Assemble the system matrix
    K = stiffness(femm, sysmat_assembler_sparse, geom, u,0*u);;
    
    % Load
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[-Inf,Inf,-Inf,Inf,Length,Length],'inflate',Thickness/1000));
    lfemm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',subset(bdry_fes,bcl),...
        'integration_rule',gauss_rule(struct('dim',2,'order',2))));
    fi= force_intensity(struct ('magn',[0;0;-magn; ]));
       F =  distrib_loads(lfemm, sysvec_assembler, geom, u, fi, 2);
    % Solve
    u = scatter_sysvec(u, K\F);
    
    
    ulin = u; %  linear displacement solution
    gv=graphic_viewer;
    cmap = cadcolors2;
    gv=reset(clear(gv,[]),[]);
    fld = field_from_integration_points(femm, geom, u, [], 'Cauchy', 3);
    nvals=fld.values; %[min(nvals),max(nvals)]
    dcm=data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
    draw(femm, gv, struct ('x',geom,'u', uscale*u, 'colorfield',colorfield, 'shrink',1));
    
%     dcm=data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
%     % Derived from element displacements
%     for i=1:count(fes)
%         femm1 = femmlock_defor_ss (struct('mater',mater,'fes',subset(fes,i),...
%             'integration_rule',gauss_rule(3, 2)));
%         fld = field_from_integration_points(femm1, geom, u, [], 'Cauchy', 3);
%         nv=get(fld,'values');
%         ix=find(nv~=0);
%         nv(ix)=mean(nv(ix));
%         colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, nv)));
%         draw(femm1, gv, struct ('x',geom,'u', uscale*u, 'colorfield',colorfield, 'shrink',1));
%     end
%         
draw_colorbar(gv, struct('colormap',cmap,'position',[0.85 0.15 0.05 0.7],...
                'minmax',dcm.range,'label',['\sigma_{z}']));
    pause(0.01); % next eigenvector
    
    
    
    disp('Stability analysis');
    [~,femm1] = restoring_force(femm,sysvec_assembler, geom,ulin,0*ulin);       % Internal forces
    K_geo = stiffness_geo(femm1, sysmat_assembler_sparse, geom, ulin, ulin);            % Geometric stiffness
    
    options.tol =1e-26;
    options.maxit = 5000;
    options.disp = 0;
    [W,Omega]=eigs(-(K_geo+K_geo')/2,(K+K')/2,...
        neigvs,'lm', options);
    Omegas=(1./(diag(Omega)));
   
    
    gv=graphic_viewer;
    cmap = cadcolors2;
    w=clone(ulin,'w');
    for i=neigvs:-1:1
        disp(['  Eigenvector ' num2str(i) ' eigenvalue ' num2str(Omegas(i)) ', error ' num2str(abs((Omegas(i)-Reference(i))/Reference(i)*100)) '%']);
        gv=reset(clear(gv,[]),[]);
        
        w = scatter_sysvec(w, W(:,i));
        v = w.values;
        wmag = nodal_field(struct ('name',['wmag'], 'dim', 1, 'nfens',count(fens)));
        wmag = scatter(wmag,1:wmag.nfens,sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2));
        nvals =gather_values(wmag,1:wmag.nfens);
        nvalsrange=[min(nvals),max(nvals)];
        dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
        draw(femm,gv, struct ('x', geom, 'u', scale*w,'colorfield',colorfield, 'shrink',1.0));
        set_graphics_defaults
        pause(1); % next eigenvector
    end
    %     assignin('caller','faesor_test_passed',(norm(Omegas(1)-26.038696877297280)<1e-3));
end
