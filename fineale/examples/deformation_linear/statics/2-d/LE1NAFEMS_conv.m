% LE1 NAFEMS benchmark, convergence analysis
disp('LE1 NAFEMS benchmark, convergence analysis');

% Parameters:
E = 210e3;% 210e3 MPa
nu = 0.3;
p = 10;% 10 MPA Outward pressure on the outside ellipse
sigma_yD= 92.7;% tensile stress at [2.0, 0.0] meters
graphics=0;

integration_rule = trapezoidal_rule(struct('dim',2));
Edge_integration_rule =trapezoidal_rule(struct('dim',1));
Convertf=[];
femmf =@femm_deformation_linear;
Style ='ks-'; Label='Q4';

integration_rule = simpson_1_3_rule(struct('dim',2));
Edge_integration_rule =simpson_1_3_rule(struct('dim',1));
Convertf=@Q4_to_Q8;
femmf =@femm_deformation_linear;
Style ='ks-'; Label='Q8';

ns=[5,10,20,40]; Results = [];
% ns=[4,8,16]; Results = [];
for n=ns;% number of elements through the thickness
    
    
    % Mesh'
    [fens,fes]=Q4_block(1.0,pi/2, n, n*2, struct('other_dimension', 0.1*1000));
    if (~isempty( Convertf ))
        [fens,fes] = Convertf(fens,fes, struct('other_dimension', 0.1*1000));
    end
    bdry_fes = mesh_boundary(fes, struct('other_dimension', 0.1*1000));
    icl = fe_select(fens, bdry_fes, struct('box', [1,1,0,pi/2],'inflate',1/100));
    for i=1:count(fens)
        t=fens.xyz(i,1); a=fens.xyz(i,2);
        fens.xyz(i,:)=1000*[(t*3.25+(1-t)*2)*cos(a) (t*2.75+(1-t)*1)*sin(a)];
    end
    % drawmesh( {fens,fes})
    % view (2)
    % return
    % Material
    prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
    mater = material_deformation_linear_biax (struct('property',prop, ...
        'reduction','stress'));
    % Finite element block
    femm = femmf(struct ('material',mater, 'fes',fes, 'integration_rule',integration_rule));
    efemm = femmf (struct ('material',mater, 'fes',subset(bdry_fes,icl),...
        'integration_rule',Edge_integration_rule));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
    % Define the displacement field
    u   = clone(geom,'u');
    u   = u*0; % zero out
    % Apply EBC's
    % plane of symmetry
    ebc_fenids=fenode_select (fens,struct('box',[0 Inf 0 0],'inflate',1/10000));
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fixed*0+2;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    % plane of symmetry
    ebc_fenids=fenode_select (fens,struct('box',[0 0 0 Inf],'inflate',1/10000));
    ebc_fixed=ones(1,length (ebc_fenids));
    ebc_comp=ebc_fixed*0+1;
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    % tic
    % Assemble the system matrix
    K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
    % Load
    fi=force_intensity(struct('magn',@(x) (p*[2.75/3.25*x(1),3.25/2.75*x(2)]'/norm([2.75/3.25*x(1),3.25/2.75*x(2)]))));
    F = distrib_loads(efemm, sysvec_assembler, geom, u, fi, 2);
    % Solve
    u = scatter_sysvec(u, K\F);
    
    % Plot
    if graphics
        gv=graphic_viewer;
        gv=reset (gv,struct ('limits',1000*[0 1.06*3.25 0 1.6*3.25]));
        scale=1000;
        cmap = jet;
    end
    fld = field_from_integration_points(femm, geom, u, [], 'Cauchy',2);
    corner=fenode_select (fens,struct('box',1000*[2.0 2.0 0 0],'inflate',1/10000));
    Approx_sigma_yD =gather_values (fld, corner);
    disp( ['n=' num2str(n) ',  Corner stress sigma_y=' num2str(Approx_sigma_yD) ' MPa'])
    disp( [' i.e. ' num2str(gather_values(fld, corner)/sigma_yD*100) '% of the Reference value'])
    Results = [Results,Approx_sigma_yD];
    if graphics
        nvals=fld.values;
        nvalsrange=[min(nvals),max(nvals)];
        dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
        draw(femm,gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink',1.0));
        draw(femm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
        % draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
        colormap(cmap);
        cbh=colorbar;
        set(cbh,...
            'Position',[0.815 0.15 0.05 0.7],...
            'YLim',[0,1],...
            'YTick',[0,1],...
            'YTickLabel',{[num2str(nvalsrange(1))],[num2str(nvalsrange(2))]});%{[num2str((min(nvals)))],[num2str((max(nvals)))]}
        set(get(cbh,'XLabel'),'String','\sigma_y');
        view (2)
    end
end

% integration_rule = simpson_1_3_rule(2);
% Edge_integration_rule =simpson_1_3_rule(1);
% Convertf=@Q4_to_Q8;
% femmf =@femmlock_defor_ss;
% Style ='ks-'; Label='Q8';
% ns= [3     6    12    24    48    96]
% Results= [  89.875976098320010  92.853550715984213  93.029698724562550  92.814638226531883  92.706647258582095] 92.671530936023160
[xestim, beta, C] = richextrapol( Results(end-2:end),1./ns(end-2:end))
loglog(1./ns,abs(Results-xestim)/sigma_yD, 'r-o')
grid on
