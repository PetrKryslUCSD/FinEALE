% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Concrete column with temperature boundary condition.  Finite element mesh
% that consists of a single four-node quadrilateral.
function pnpConcreteColumn1Q4

%% 
% Evaluate for a given data
a=2.5; dy=a/2*sin(15/180*pi); dx=a/2*cos(15/180*pi); Q=4.5; k=1.8; Dz=1.0;
x=[0,0; 2*dx,-2*dy; a,0; 2*dx,2*dy];  % Coordinates
%% 
% 
%Define a little function to calculate the gradient of the basis functions
%in the parametric coordinates.
gradNpar =@(xi,eta)[ 
[   eta/4 - 1/4,   xi/4 - 1/4]
[   1/4 - eta/4, - xi/4 - 1/4]
[   eta/4 + 1/4,   xi/4 + 1/4]
[ - eta/4 - 1/4,   1/4 - xi/4]];
%The expressions actually come from  the following symbolic code:
% %Define the basis function matrix
% syms  xi eta real
% N=[(xi-1)*(eta-1)/4;(xi+1)*(eta-1)/-4;(xi+1)*(eta+1)/4;(xi-1)*(eta+1)/-4]
% % Differentiate to obtain the Basis function gradients
% gradNpar=[diff(N,'xi'),diff(N,'eta')]

%% 
% These are the integration point data
xe=[-0.577350269189626,-0.577350269189626;
    -0.577350269189626,+0.577350269189626;
    +0.577350269189626,-0.577350269189626;
    +0.577350269189626,+0.577350269189626];
W=[1,1,1,1];

%% 
% Initialize the elementwise conductivity matrix.
K=zeros(4);
%% 
% Loop over the quadrature points.
for q=1:4
    xi=xe(q,1);     eta=xe(q,2);
    %Compute the Jacobian matrix
    J=(x'*gradNpar(xi,eta));
    % The Jacobian
    detJ=det(J)
    % The  gradient  of the basis functions with respect to x,y 
    gradN= gradNpar(xi,eta)/J;
    % Add the contribution to the conductivity matrix
    K=K+k*Dz*gradN*gradN'*detJ*W(q);
end

matrix2latex(K,struct('precision' ,4, 'separator',','))

%% 
% We will find it convenient to define a little function to evaluate the
% basis function values at a given quadrature point location.
N=@(xi,eta)[(xi-1)*(eta-1)/4;(xi+1)*(eta-1)/-4;(xi+1)*(eta+1)/4;(xi-1)*(eta+1)/-4];

%% 
% Initialize the elementwise heat-load vector .
L=zeros(4,1);
%% 
% Loop over the quadrature points.
for q=1:4
    xi=xe(q,1);     eta=xe(q,2);
    %Compute the Jacobian matrix
    J=(x'*gradNpar(xi,eta));
    % The Jacobian
    detJ=det(J);
    % Add the contribution to the conductivity matrix
    L=L+Q*Dz*N(xi,eta)*detJ*W(q);
end

%% 
% The solution vector consists of all zeros, except in the first entry.
T=zeros(4,1);
%Solve the global equations
T(1)=K(1,1)\L(1)
%% 
% Now we will post process to extract the heat flux vectors at the
% quadrature points.
%% 
% Loop over the quadrature points.
for q=1:4
    xi=xe(q,1);     eta=xe(q,2);
    % Compute the location of the quadrature point
    qp=N(xi,eta)'*x
    %Compute the Jacobian matrix
    J=(x'*gradNpar(xi,eta));
    % The  gradient  of the basis functions with respect to x,y 
    gradN= gradNpar(xi,eta)/J;
    % Compute the gradient of temperature
    gradT  =T'*gradN
    % Heat flux vector
    q=-k*gradT'
end

end