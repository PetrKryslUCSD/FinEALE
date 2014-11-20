function sw129b_l3_irreg
% Taut wire: statics,  triangular load (128-b);
L=6;
P=4;
q0= 0.1;
q =@(x)(q0*(L-x)/L);% This is consistent now with the sign convention in the textbook
n=3; % number of elements
tolerance =L/n/100;
    % Mesh
[fens,fes]= L2_block(L,n, 1.0);
nen=count(fens);
[fens,fes] = L2_to_L3(fens,fes, []);
x=fens.xyz;
x(nen+1:end,:)=x(nen+1:end,:)+(-1)*L/n/6;
fens.xyz=x;
% Finite element block
prop = property_deformation_linear_iso (struct ('E',P,'nu',0,'rho',1.0));
    mater = material_deformation_linear_uniax (struct('property',prop));
    feb = femm_deformation_linear(struct ('material',mater,...
        'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',1,'order',4))));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 1, 'fens',fens));
% Define the displacement field
w   = 0*clone(geom,'w');
% Apply EBC's
n0= fenode_select(fens,struct('distance',tolerance, 'from',0));
nL= fenode_select(fens,struct('distance',tolerance, 'from',L));
fenids=[n0,nL]; prescribed=[1,1]; component=[1,1]; val=[0,0];
w   = set_ebc(w, fenids, prescribed, component, val);
w   = apply_ebc (w);
% Number equations
w   = numberdofs (w);
% Assemble the system matrix
K = stiffness(feb, sysmat_assembler_sparse, geom, w);
% Load
fi = force_intensity(struct ('magn',q));
F = distrib_loads(feb, sysvec_assembler, geom, w, fi,3);
% Solve
w = scatter_sysvec(w, K\F);
% Plot

figure; hold on;
plot(geom.values,0*geom.values,'ko-','linewidth',3);
xlabel('x'),ylabel('Mesh')

figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    ix=isoparamipol1d(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))));
    iw=isoparamipol1d(fes,pc,gather_values(w,connected_nodes(subset(fes,j))));
    plot(ix/L,iw,'r-','linewidth',3);hold on
    x=ix;
    plot(x/L,(q0*x.^3)/(6*L*P) - (q0*x.^2)/(2*P) + (L*q0*x)/(3*P),'k:','linewidth',3);
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('w')

figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    ix=isoparamipol1d(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))));
    iw=isoparamipol1d(fes,pc,gather_values(w,connected_nodes(subset(fes,j))));
    x=ix;
    plot(ix/L,(q0*x.^3)/(6*L*P) - (q0*x.^2)/(2*P) + (L*q0*x)/(3*P)-iw,'r-','linewidth',3);hold on
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('w_{ex}-w')


figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    ix=isoparamipol1d(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))));
    iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))),gather_values(w,connected_nodes(subset(fes,j))));
    plot(ix/L,P*iwp,'r-','linewidth',3);hold on
    x=ix;
    plot(x/L,P*((L*q0)/(3*P) - (q0*x)/P + (q0*x.^2)/(2*L*P)),'k:','linewidth',3);
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('Pw^{\prime}')

figure; hold on;
pc=linspace(-1,+1,100);
for j=1:count(fes)
    ix=isoparamipol1d(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))));
    iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))),gather_values(w,connected_nodes(subset(fes,j))));
    x=ix;
    plot(ix/L,P*((L*q0)/(3*P) - (q0*x)/P + (q0*x.^2)/(2*L*P))-P*iwp,'r-','linewidth',3);hold on
end
left_handed_axes;% Plotting is consistent with the textbook sign convention
xlabel('x/L'),ylabel('Pw_{ex}^{\prime}-Pw^{\prime}')

disp('Reaction at x=0')
q(0)*L/2*(2*L/3)/L

disp(' Finite element result for the reaction')
j=1;
pc=-1;
ix=isoparamipol1d(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))));
    iwp=isoparamipol1ddersp(fes,pc,gather_values(geom,connected_nodes(subset(fes,j))),gather_values(w,connected_nodes(subset(fes,j))));
reaction=P*iwp
end

