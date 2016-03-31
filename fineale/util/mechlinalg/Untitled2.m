% U=diag(rand(3,1));
% U = (U + U')/2;      % Force Hermitian by taking nearest Hermitian matrix.
% [R,A]=qr(rand(3)); det(R)
% F=R*U;
% [R1,U1] = polardecomp(F)
% U-U1
% R-R1

function f
A=rand(3);
B=rand(3);
tensor_3x3t_double_contraction(A,B)
v=0;
for    i =1:3
for    j =1:3
v=v+A(i,j)*B(i,j);
end; clear j
end; clear i
v
function v = tensor_3x3t_double_contraction(A,B)
            v = sum(sum(A.*B));
end
end