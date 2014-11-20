% Solution of the same problem as in pnpSquareInSquare but with ad hoc
% Matlab code.
function pnpSingleConductivity1
% Coordinates of nodes
x= [-1,0;  1/2,-(3^(1/2))/2; 1/2,(3^(1/2))/2];
% Connectivity
conn = [1,2,3];
Dz=1; % Thickness of the slice assumed equal to 1.0!
kappa=[2.0 0; 0 0.1]; % orthotropic conductivity matrix
alpha =90;% local material orientation angle in degrees



Nder= [-1,-1;1,0;0,1];% gradients of the basis functions WRT parametric
J=x(conn(1,:),:)'*Nder;% compute the Jacobian matrix
Ndersp=Nder/J;%  compute the spatial gradients of the basis functions
% Elementwise conductivity matrix: note  the material orientation
% matrix
Rm = [cos(2*pi/360*alpha), -sin(2*pi/360*alpha); ...
      sin(2*pi/360*alpha), cos(2*pi/360*alpha)];% local material directions
Ke=  double(det(J)/2*Ndersp*Rm*kappa*Rm'*Ndersp')*Dz

% Heat load vector due to unit drop in temperature in the X direction
Ke*[1;0;0]

% Draw the temperature distribution  on the mesh
patch('xdata',x(conn(1,:),1),'ydata',x(conn(1,:),2), 'facecolor','r');;
for j=1:size(x,1)
    text(x(j,1),x(j,2),num2str(j),'backgroundcolor','w', 'fontname', 'times', 'fontsize',18)
end 
axis equal
xlabel x, ylabel y, zlabel T 
set(gca,'fontname', 'times','fontsize',18)
set(get(gca,'xlabel'),'fontname', 'times','fontsize',18)
set(get(gca,'ylabel'),'fontname', 'times','fontsize',18)
grid on