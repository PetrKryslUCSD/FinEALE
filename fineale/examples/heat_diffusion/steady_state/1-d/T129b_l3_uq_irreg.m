% Wall with uniform distributed source and quadratic elements.
function T129b_l3_uq_irreg
graphics =~0;
L=6;
P=4;
q0= 0.1;
q =@(x)(q0);% This is consistent now with the sign convention in the textbook
n=3; % number of elements
node_shift =(+1)*L/n/5;
% Mesh
[fens,fes]= L2_block(L,n, 1.0);
nen=count(fens);
[fens,fes]= L2_to_L3(fens,fes, 1.0);
x=fens.xyz;
x(nen+1:end,:)=x(nen+1:end,:)+node_shift;
fens.xyz=x;
% Finite element block
prop=property_heat_diffusion(struct('thermal_conductivity',P,'source',0));
mater=material_heat_diffusion (struct('property',prop));
femm = femm_heat_diffusion (struct (...
    'material',mater,'fes',fes,...
    'integration_rule',gauss_rule(struct( 'dim',1,'order',3))));
geom = nodal_field(struct ('name',['geom'], 'dim', 1, 'fens',fens));
tempn = nodal_field(struct ('name',['temp'], 'dim', 1,...
    'nfens',geom.nfens));
Fenids=[fenode_select(fens,struct ('box',[0,0],  'inflate',L/1000)),fenode_select(fens,struct ('box',[L,L],  'inflate',L/1000))];
tempn = set_ebc(tempn, Fenids, Fenids*0+1, Fenids*0+1, Fenids*0);
tempn = apply_ebc (tempn);
tempn = numberdofs (tempn);
 

% Assemble the system matrix
K = conductivity(femm, sysmat_assembler_sparse, geom, tempn);
% Load
fi= force_intensity(struct('magn',q));
F = distrib_loads(femm, sysvec_assembler, geom, tempn, fi, 3);
% Solve
w = scatter_sysvec(tempn, K\F);
% Plot
if (graphics)
    
figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    fe=subset(fes,j);
    ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
    iw=isoparamipol1d(fes,pc,gather_values(w,fe.conn));
    plot(ix/L,iw,'r-','linewidth',3);hold on
end
x=(pc+1)/2*L;
plot(x/L,q(0)/(2*P)*x.*(L-x),'k-','linewidth',3);hold on
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('w')

figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
   fe=subset(fes,j);
    ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
    iw=isoparamipol1d(fes,pc,gather_values(w,fe.conn));
    x=ix;
    plot(ix/L,q(0)/(2*P)*x.*(L-x)-iw,'r-','linewidth',3);hold on
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('w_{ex}-w')


figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    fe=subset(fes,j);
    ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
    iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,fe.conn),gather_values(w,fe.conn));
    plot(ix/L,P*iwp,'r-','linewidth',3);hold on
    x=ix;
    plot(x/L,P*(q(0)/(2*P)*(L-2*x)),'k:','linewidth',3);
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('Pw^{\prime}')

figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    fe=subset(fes,j);
    ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
    iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,fe.conn),gather_values(w,fe.conn));
    x=ix;
    plot(ix/L,P*(q(0)/(2*P)*(L-2*x))-P*iwp,'r-','linewidth',3);hold on
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('Pw_{ex}^{\prime}-Pw^{\prime}')

disp('Reaction at x=0')
q(0)*L/2*(2*L/3)/L

disp(' Finite element result for the reaction')
j=1;
pc=-1;
fe=subset(fes,j);
    ix=isoparamipol1d(fes,pc,gather_values(geom,fe.conn));
iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,fe.conn),gather_values(w,fe.conn));
P*iwp

end

