% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Concrete column with hydration heat and convection condition on the wet boundary.
function pnpConcreteColumnsymb_conv

%%
% Define some useful quantities.  Some of the calculations are going to be
% done in symbolic variables.
syms k dx dy Q Dz real % the variables in the problem, k is k
a=2.5; dy=a/2*sin(15/180*pi); dx=a/2*cos(15/180*pi); Q=4.5; k=1.8; Dz=1.0;
h= 5.;
gradNpar= [-1,-1;1,0;0,1];%  Gradients of the basis functions wrt the parametric coords
xall= [0,0; dx,-dy; dx,dy; 2*dx,-2*dy; 2*dx,2*dy];%Coordinates of the nodes
dof=[3    1    2    5    4];% Numbers of the degrees of freedom
N_f=5;
%%
% Global conductivity matrix and heat load vector.
K=sym(zeros(5));
L=sym(zeros(5,1));
%%
%
%First element
conn= [1,2,3];% The definition of the element, listing its nodes
x=xall(conn,:);% The coordinates  of the three nodes
J=x'*gradNpar % Compute the Jacobian matrix
Se=det(J )/2 % The area of the triangle
gradN=gradNpar/J
(Se*gradN(1,:)*gradN(1,:)'*k*Dz)
(Se*gradN(1,:)*gradN(2,:)'*k*Dz)
Ke1= (Se*gradN*gradN'*k*Dz)
edof=dof(conn)
K(edof,edof)=K(edof,edof)+Ke1;
LeQ1=Se*Q*Dz/3*ones(3,1);
L(edof)=L(edof)+LeQ1;



%%
%
%Second element
conn= [2,4,5];
x=xall(conn,:);
J=x'*gradNpar
Se=det(J )/2
gradN=gradNpar/J
Ke2= (Se*gradN*gradN'*k*Dz)
edof=dof(conn)
K(edof,edof)=K(edof,edof)+Ke2;
LeQ2=Se*Q*Dz/3*ones(3,1);
L(edof)=L(edof)+LeQ2;

%%
%
%Third element
conn= [2,5,3];
x=xall(conn,:);
J=x'*gradNpar
Se=det(J )/2
gradN=gradNpar/J
edof=dof(conn)
Ke3= (Se*gradN*gradN'*k*Dz)
K(edof,edof)=K(edof,edof)+Ke3;
LeQ3=Se*Q*Dz/3*ones(3,1);
L(edof)=L(edof)+LeQ3;


%% 
% The surface heat transfer matrix for element [4,5]
conn= [4,5];
he=norm(diff(xall(conn,:)))
He=h*Dz*he/6*[2,1;
              1,2]
%Surface heat transfer matrix with the trapezoidal rule is
%diagonal
% He=h*Dz*he/2*[1,0;
%               0,1]
edof=dof(conn)      
K(edof,edof)=K(edof,edof)+He;
%% 
% Compute the solution.
T=eval(K)\L;

%% 
% Substitute numerical values and evaluate the result.
m=vpa(eval(T),4)
% matrix2latex(double(m),struct('precision' ,5, 'separator',','))
