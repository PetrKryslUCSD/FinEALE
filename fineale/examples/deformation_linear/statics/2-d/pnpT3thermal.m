% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Thermal loading nodal forces for a single triangle.
function pnpT3thermal

%% 
xall= [-1, -1/2; 3, 2; 1, 2];%Coordinates of the nodes

gradNpar= [-1,-1;1,0;0,1];%Gradients of the basis fncs wrt the param. coords
conn= [1,2,3];% The definition of the element, listing its nodes
x=xall(conn,:);% The coordinates  of the three nodes
J=x'*gradNpar % Compute the Jacobian matrix
gradN=gradNpar/J

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

syms DeltaT CTE E nu t real 
D=E/(1-nu^2)*[1,nu,0; nu,1,0; 0,0,(1-nu)/2]
eth=DeltaT*CTE*[1;1;0];
Se=det(J)/2
F1=simplify(Se*t*B1'*D*eth)
F2=simplify(Se*t*B2'*D*eth)
F3=simplify(Se*t*B3'*D*eth)

figure
patch('xdata',x(conn,1),'ydata',x(conn,2), 'facecolor','y', 'linewidth',3);;
grid on
axis equal
set(gca,'xlim',[-2,5])
set(gca,'ylim',[-2,5])
labels('\it x','\it y')
% ah=annotation('arrow',[.9 .5],[.9,.5],'Color','r');
line([x(1,1),x(1,1)+F1(1)/(CTE*DeltaT*E*t/(1-nu))],[x(1,2),x(1,2)+F1(2)/(CTE*DeltaT*E*t/(1-nu))],'Color','r', 'linewidth',4)
line([x(2,1),x(2,1)+F2(1)/(CTE*DeltaT*E*t/(1-nu))],[x(2,2),x(2,2)+F2(2)/(CTE*DeltaT*E*t/(1-nu))],'Color','r', 'linewidth',4)
line([x(3,1),x(3,1)+F3(1)/(CTE*DeltaT*E*t/(1-nu))],[x(3,2),x(3,2)+F3(2)/(CTE*DeltaT*E*t/(1-nu))],'Color','r', 'linewidth',4)
set_pub_defaults