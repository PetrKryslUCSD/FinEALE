% Square domain with internal heat source  and fixed temperature 
% on the boundary. Ad hoc symbolic solution.
function pnpSquare1symb
syms a real
Nder= [-1,-1;1,0;0,1];
x= [0,0;a,0;0,a];
J=x'*Nder
det(J )
Ndersp=Nder/J
syms Dz kappa Q real
K3= (a^2/2*Ndersp*Ndersp'*kappa*Dz)
K=sym(zeros(6 ));
eq= [4,5,1];
K(eq,eq)=K(eq,eq)+K3;
eq= [5,6,2];
K(eq,eq)=K(eq,eq)+K3
eq= [1,2,3];
K(eq,eq)=K(eq,eq)+K3
eq= [2,1,5];
K(eq,eq)=K(eq,eq)+K3
matrix2latex(K/(Dz*kappa))

K=K(1:3,1:3)
L=-20*Dz*kappa*[-1/2;-1;0]
K\L

L=Q*(a^2)/2*Dz*[1;1;1/3]
K\L