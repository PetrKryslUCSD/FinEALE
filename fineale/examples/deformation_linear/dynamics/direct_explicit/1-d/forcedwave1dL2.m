function forcedwave1dL2
    % Harmonically forced wave along a prestressed cable.
    %   Detailed explanation goes here
    npart= 59;
    E=200000; rho=7.85e-9; A=200; L=7850; 
    omega= 10000;
    tolerance =L/npart/100;
    cd=sqrt(E/rho);
    dt=L/npart/cd
    nsteps = round(12.3*L/(cd)/dt);
    tspan= [0,nsteps*dt];
    tfinal =tspan(2);
    
    graphics = ~false;
    scale =1;
    
    [fens,fes]= L2_block(L,npart,A);
    n0= fenode_select(fens,struct('distance',tolerance, 'from',0));
    % Finite element block
    prop = property_deformation_linear_iso (struct ('E',E,'nu',0,'rho',rho));
    mater = material_deformation_linear_uniax (struct('property',prop));
    feb = femm_deformation_linear(struct ('material',mater,...
        'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',1,'order', 2))));
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 1, 'fens',fens));
    % Define the displacement field
    u   = 0*clone(geom,'u');
    % Apply EBC's
    % fenids=[1,n+1]; prescribed=[1, 1]; component=[1, 1]; val=[0,0];
    % u   = set_ebc(u, fenids, prescribed, component, val);
    % u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    v = 0*u;% zero velocity
    % Assemble the stiffness matrix
    K = stiffness(feb, sysmat_assembler_sparse, geom, u);
    % Assemble the consistent or lumped mass matrix
    M = lumped_mass(feb, sysmat_assembler_sparse, geom, u);
    % Solve
    U0 = gather_sysvec(u);
    V0 = gather_sysvec(v);
    A0 =M\(-K*U0);
    
    nd0 =gather_dofnums(u,n0);
    F =0*U0; F(nd0) =1000;
    
    t=0;
    o2=eigs(K,M,1,'LM');
    dt= 0.999* 2/sqrt(o2)
    
    ngraphics=  max (round(5*(1e-5)/dt),1);
    igraphics =ngraphics-1;
    if graphics
        set_graphics_defaults(figure)
    end
    u0s=[];
    while t <tfinal
        U1 = U0 +dt*V0+(dt^2)/2*A0;
        A1 = M\(-K*U1 +F*sin(omega*t));
        V1 = V0 +dt/2* (A0+A1);
        U0 = U1;
        V0 = V1;
        A0 = A1;
        t=t+dt;
        igraphics=igraphics+1;
        if graphics & igraphics==ngraphics
            u = scatter_sysvec(u, U1);
            X =gather_values(geom,1:count(fens));
            U =gather_values(u,1:count(fens));
            plot(X,U);
            labels('Distance along cable x','Deflection u')
            set(gca,'ylim',[-0.,3*0.6]);
            set(gca,'xlim',[0,L]);
            view(2); pause(0.1);
            igraphics= 0;
        end
        u0s=[u0s,U1(nd0)];
    end
    
    figure;
    plot(u0s);
end

