% Derivation of the transformation matrix [T]
clear all;
syms T alpha R
syms a a11 a12 a13 a21 a22 a23 a31 a32 a33
a = [a11,a12,a13;
     a21,a22,a23;
     a31,a32,a33];
 % % it can be done in terms of l,m,n's as well
 % syms a l1 m1 n1 l2 m2 n2 l3 m3 n3
 % a = [l1,m1,n1;l2,m2,n2;l3,m3,n3]
T(1:6,1:6) = 0;
for i=1:1:3
for j=1:1:3
 if i==j; alpha = j; else alpha = 9-i-j; end
 for p=1:1:3
 for q=1:1:3
  if p==q beta = p; else beta = 9-p-q; end
  T(alpha,beta) = 0;
  if alpha<=3 & beta<= 3; T(alpha,beta)=a(i,p)*a(i,p); end
  if alpha> 3 & beta<= 3; T(alpha,beta)=a(i,p)*a(j,p); end
  if alpha<=3 & beta>3; T(alpha,beta)=a(i,q)*a(i,p)+a(i,p)*a(i,q);end
  if alpha>3 & beta>3; T(alpha,beta)=a(i,p)*a(j,q)+a(i,q)*a(j,p);end
 end
 end
end
end
T
R = eye(6,6); R(4,4)=2; R(5,5)=2; R(6,6)=2; % Reuter matrix
Tbar = R*T*R^(-1)