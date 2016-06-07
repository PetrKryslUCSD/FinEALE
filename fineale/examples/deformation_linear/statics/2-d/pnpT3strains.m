% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Strain patterns calculated for a single triangle.
function pnpT3strains

%% 
xall= [-1, -1/2; 3, 2; 1, 2];%Coordinates of the nodes

gradNpar= [-1,-1;1,0;0,1];%Gradients of the basis fncs wrt the param. coords
conn= [1,2,3];% The definition of the element, listing its nodes
x=xall(conn,:);% The coordinates  of the three nodes
J=x'*gradNpar % Compute the Jacobian matrix
gradN=gradNpar/J

figure
patch('xdata',x(conn,1),'ydata',x(conn,2), 'facecolor','y', 'linewidth',3);;
grid on
axis equal
set(gca,'xlim',[-1,3])
set(gca,'ylim',[-1,2])
labels('\it x','\it y')
set_pub_defaults
%% 
% Anonymous function to calculate one nodal strain-displacement matrix
Bn=@(gradNn)[gradNn(1),0;
    0,gradNn(2);
    gradNn(2),gradNn(1)];

%% 
% 
B1=Bn(gradN(1,:))
B2=Bn(gradN(2,:))
B3=Bn(gradN(3,:))

syms psi Xc Yc real 
syms x1 y1 x2 y2 x3 y3 real 
rot=@(X)[0,-psi;psi,0]*(X-[Xc;Yc]);
u1=rot([x1;y1])
u2=rot([x2;y2])
u3=rot([x3;y3])
e=B1*u1+B2*u2+B3*u3
simplify(e)