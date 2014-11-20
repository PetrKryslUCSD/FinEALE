% Wall with uniform distributed load (Exercise 129-b).
function T129b_l2_uq
graphics =~0;
L=6;
kappa=4;
Q =.1;% 
n=7; % number of elements

[fens,fes]= L2_block(L,n, 1.0);% Mesh
% Finite element block
prop=property_heat_diffusion(struct('thermal_conductivity',kappa,'source',Q));
mater=material_heat_diffusion (struct('property',prop));
femm = femm_heat_diffusion (struct (...
    'material',mater,'fes',fes,...
    'integration_rule',gauss_rule(struct( 'dim',1,'order',2))));
geom = nodal_field(struct ('name',['geom'], 'dim', 1, 'fens',fens));
tempn = nodal_field(struct ('name',['temp'], 'dim', 1,...
    'nfens',geom.nfens));
Fenids=[1,count(fens)];
tempn = set_ebc(tempn, Fenids, Fenids*0+1, Fenids*0+1, Fenids*0);
tempn = apply_ebc (tempn);
tempn = numberdofs (tempn);
 

% Assemble the system matrix
K = conductivity(femm, sysmat_assembler_sparse, geom, tempn);
% Load
fi= force_intensity(struct('magn',Q));
F = distrib_loads(femm, sysvec_assembler, geom, tempn, fi, 3);
% Solve
tempn = scatter_sysvec(tempn, K\F);
% Plot
if (graphics)
    figure; hold on;
    pc=linspace(-1,+1,100);
    for j=1:count(fes)
        fe=subset(fes,j);
        ix=isoparamipol1d(fe,pc,gather_values(geom,fe.conn));
        iw=isoparamipol1d(fe,pc,gather_values(tempn,fe.conn));
        plot(ix/L,(kappa/Q/L^2)*iw,'r-','linewidth',3);hold on
    end
    x=(pc+1)/2*L;
    plot(x/L,(kappa/Q/L^2)*Q/(2*kappa)*x.*(L-x),'k:','linewidth',3);hold on
    xlabel('x/L'),ylabel('(\kappa/Q/L^2)*T')
    
    figure; hold on;
    pc=linspace(-1,+1,100);
    for j=1:count(fes)
        fe=subset(fes,j);
        ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
        iw=isoparamipol1d(fes,pc,gather_values(tempn,fe.conn));
        x=ix;
        plot(ix/L,(kappa/Q/L^2)*Q/(2*kappa)*x.*(L-x)-(kappa/Q/L^2)*iw,'r-','linewidth',3);hold on
    end
    xlabel('x/L'),ylabel('(\kappa/Q/L^2)*(T_{ex}-T)')
    
    figure; hold on;
    pc=linspace(-1,+1,10);
    for j=1:count(fes)
        fe=subset(fes,j);
        ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
        iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,fe.conn),gather_values(tempn,fe.conn));
        plot(ix/L,-kappa*iwp,'r-','linewidth',3);hold on
    end
    x=(pc+1)/2*L;
    plot(x/L,-kappa*Q/2/kappa*(L-2*x),'k:','linewidth',3);hold on
    xlabel('x/L'),ylabel('q')

    figure; hold on;
    pc=linspace(-1,+1,10);
    for j=1:count(fes)
        fe=subset(fes,j);
        ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
        iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,fe.conn),gather_values(tempn,fe.conn));
        x=ix;
        plot(ix/L,-kappa*Q/2/kappa*(L-2*x)-(-kappa*iwp),'r-','linewidth',3);hold on
    end
%     left_handed_axes;% Plotting is consistent with the textbook sign convention
    xlabel('x/L'),ylabel('q_{ex}-q')
end

