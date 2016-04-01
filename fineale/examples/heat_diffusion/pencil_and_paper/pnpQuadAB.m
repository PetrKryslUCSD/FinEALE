function pnpQuadAB

syms A B real
x=[A,B;0,B;0,0;A,0];
%Define the basis function matrix
syms  xi eta real
N=[(xi-1)*(eta-1)/4;(xi+1)*(eta-1)/-4;(xi+1)*(eta+1)/4;(xi-1)*(eta+1)/-4]
% Differentiate to obtain the Basis function gradients
gradNpar=[diff(N,'xi'),diff(N,'eta')]
%Compute the Jacobian matrix
J=simplify(x'*gradNpar)
% The Jacobian
det(J)
% The  gradient  of the basis functions with respect to x,y is a
% symbolic expression which needs to be evaluated for particular values of
% xi,eta
gradN= gradNpar/J
%Note that using subs() will substitute the values of the parametric
%coordinates
xi=-0.577350269189626; eta=-0.577350269189626;
K1 =subs(gradN)*subs(gradN)';
xi=-0.577350269189626; eta= 0.577350269189626;
K2 =subs(gradN)*subs(gradN)';
xi= 0.577350269189626; eta=-0.577350269189626;
K3 =subs(gradN)*subs(gradN)';
xi=0.577350269189626 ; eta=0.577350269189626;
K4 =subs(gradN)*subs(gradN)';

%Intermediate result: we still need to multiply with  the thermal
%conductivity  and the thickness of the slice, but we do include the
%constant Jacobian here
K=simplify(det(J)*(K1+K2+K3+K4))

matrix2latex((6*A*B)*K,struct('precision' ,4, 'separator',','))
end