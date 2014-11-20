% Infinite plane with a circular hole.
%
% Infinite plane with a circular hole of radius a, subject to  
% uniform tension in the X direction. Analytical solution is available.
%
% This version uses the mechanism of bundled algorithms. The solution
% is obtained with the algorithm deformation_linear_statics().
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
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.fes= fes;
    region.E =E;
    region.nu=nu;
    region.integration_rule =gauss_rule(struct('dim', 3,'order', 2));
    model_data.region{1} =region;
    
    clear essential
    essential.component= [1];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',tol));
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.component= [2];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[-Inf,Inf,0,0,-Inf,Inf],'inflate',tol));
    model_data.boundary_conditions.essential{2} = essential;
    
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0],'inflate',tol)),fenode_select(fens,struct('box',[-Inf,Inf,-Inf,Inf,H,H],'inflate',tol))];
    model_data.boundary_conditions.essential{3} = essential;
    
    bdry_fes = mesh_boundary(fes, []);
    bclx = fe_select(fens, bdry_fes, ...
        struct ('box',[L,L,-Inf,Inf,-Inf,Inf],'inflate',tol));
    clear traction
    traction.fes =subset(bdry_fes,bclx);;
    traction.traction= @(x) ([sigmaxx(x),sigmaxy(x),0]);
    traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));
    model_data.boundary_conditions.traction{1} = traction;
    
    bcly = fe_select(fens, bdry_fes, ...
        struct ('box',[-Inf,Inf,W,W,-Inf,Inf],'inflate',tol));
    clear traction
    traction.fes =subset(bdry_fes,bcly);;
    traction.traction= @(x) ([sigmaxy(x),sigmayy(x),0]);
    traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));
    model_data.boundary_conditions.traction{2} = traction;
    
    % Solve
    tic; model_data =deformation_linear_statics(model_data);
    toc
    
    
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
    
    % Plot
    gv=graphic_viewer;
    gv=reset (gv,struct ([]));
    scale=10000;
    cmap = jet;
    fld = field_from_integration_points(model_data.region{1}.femm, model_data.geom, model_data.u, [], 'Cauchy',1);
    nvals=fld.values/(U.MEGA*U.PA);%min(nvals),max(nvals)
    nvalsrange=[min(nvals),max(nvals)];
    dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
    colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
    draw(model_data.region{1}.femm,gv, struct ('x', model_data.geom, 'u', +scale*model_data.u,'colorfield',colorfield, 'shrink',1.0));
    draw(model_data.region{1}.femm,gv, struct ('x', model_data.geom, 'u', 0*model_data.u,'facecolor','none', 'shrink',1));
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