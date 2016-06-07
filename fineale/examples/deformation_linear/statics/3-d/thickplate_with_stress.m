function doit
    % Thick plate, small displacement, small strain analysis
    disp('Thick plate, small displacement, small strain analysis');

    % Parameters:
    E=1e3;
    nu=0.2;
    graphics =1; % graphic output
    uscale = 4;
    comp =6;

    % Mesh
    n=2;
    tic; % for n= 10, that is almost 14,000 equations, we need about half a minute
    [fens,fes] = H20_block(4,4,0.5,2*n,2*n,n);
    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    % Finite element block
    feb = femm_deformation_linear (struct('material',mater,'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = 0*geom; % zero out
    % Apply EBC's
    ebc_fenids=fenode_select (fens,struct('box',[-10 10 4 4 -10 10]));
    ebc_prescribed=ebc_fenids*0+1;
    ebc_comp=[];
    ebc_val=ebc_fenids*0;
    u   = set_ebc(u, ebc_fenids, ebc_prescribed, ebc_comp, ebc_val);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    % Assemble the system matrix
    K = stiffness(feb, sysmat_assembler_sparse, geom, u);
    % Load
    fi= force_intensity(struct ('magn',@(x)(-0.25*(x(1)-2)*[0, 0, 1])));
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[-10 10 -10 10 0.5 0.5],'inflate',0.001));
    lfeb = femm_deformation_linear(struct ('material',mater, 'fes',subset(bdry_fes,bcl),...
        'integration_rule',gauss_rule(struct('dim',2, 'order',2))));
    F = distrib_loads(lfeb, sysvec_assembler, geom, u, fi, 2);
    % Solve
    u = scatter_sysvec(u, K\F);
    toc

    if graphics

        figure; gv=graphic_viewer;
        gv=reset (gv,[]);
        % Continuous stress
        fld = field_from_integration_points(feb, geom, u, 0*u, 0, [], 'Cauchy', comp);
        nvals=fld.values; [min(nvals),max(nvals)]
        dcm=data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
        draw(feb, gv, struct ('x',geom,'u', uscale*u, 'colorfield',colorfield, 'shrink',1));
        camset(gv, [ 37.2735  -21.7348 25.3973    1.8154 1.7344   -0.1523  -0.4295    0.2843 0.8572    4.1432])

        figure; gv=graphic_viewer;
        gv=reset (gv,[]);
        % Ellipsoid representation of stress
        draw(feb, gv, struct ('x',geom,'u', uscale*u, 'facecolor','none'));
        id.comp=comp;
        id.container=-Inf;
        id=inspect_integration_points(feb, geom, u, 0*u, 0, [],...
            (1:count(fes)), struct ('output',['Cauchy']),...
            @mx,id)
        max_stress=id.container
        id.container=Inf;
        id=inspect_integration_points(feb, geom, u, 0*u, 0, [], ...
            (1:count(fes)), struct ('output',['Cauchy']),...
            @mn,id)
        min_stress =id.container
        dcm=data_colormap(struct ('range',[min(nvals),max(nvals)], 'colormap',jet));
        draw_integration_points(feb,gv,struct ('x',geom,'un1',u,'un',0*u,'dt',[],'u_scale',uscale,'dT',[],...
            'scale',0.05,'component',comp,'data_cmap', dcm,'cheap_arrow',0, 'linewidth',2));
        camset(gv, [ 37.2735  -21.7348 25.3973    1.8154 1.7344   -0.1523  -0.4295    0.2843 0.8572    4.1432])

        figure; gv=graphic_viewer;
        gv=reset (gv,[]);
            dcm=data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
        % Derived from element displacements
        for i=1:count(fes)
            feb = femm_deformation_linear (struct('material',mater,'fes',subset(fes,i),...
                'integration_rule',tensprod_nq_rule(struct('dim',3,'order',2))));
            fld = field_from_integration_points(feb, geom, u, 0*u, 0, [], 'Cauchy', comp);
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, fld.values)));
            draw(feb, gv, struct ('x',geom,'u', uscale*u, 'colorfield',colorfield, 'shrink',1));
        end
        camset(gv, [ 37.2735  -21.7348 25.3973    1.8154 1.7344   -0.1523  -0.4295    0.2843 0.8572    4.1432])

        figure; gv=graphic_viewer;
        gv=reset (gv,[]);
            dcm=data_colormap(struct ('range', [min(nvals),max(nvals)], 'colormap',jet));
        % Derived from element displacements
        for i=1:count(fes)
            feb = femm_deformation_linear(struct('material',mater,'fes',subset(fes,i),...
                'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
            fld = field_from_integration_points(feb, geom, u, 0*u, 0, [], 'Cauchy', comp);
            nv=fld.values;
            ix=find(nv~=0);
            nv(ix)=mean(nv(ix));
            colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nv)));
            draw(feb, gv, struct ('x',geom,'u', uscale*u, 'colorfield',colorfield, 'shrink',1));
        end
        camset(gv, [ 37.2735  -21.7348 25.3973    1.8154 1.7344   -0.1523  -0.4295    0.2843 0.8572    4.1432])

    end
end

function id= mn(id,out,xyz,u,pc)
    id.container=min(out(id.comp), id.container);
end

function id= mx(id,out,xyz,u,pc)
    id.container=max(out(id.comp), id.container);
end

