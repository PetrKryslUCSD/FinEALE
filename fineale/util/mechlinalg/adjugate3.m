function B = adjugate3(A)
% Compute the adjugate matrix of a 3 x 3 real matrix.
%
% function B = adjugate3(A)
%

% Test:
% A = rand(3)
B =0*A;
B(1,1) =(A(5)*A(9)-A(8)*A(6));
B(1,2) =-(A(4)*A(9)-A(7)*A(6));
B(1,3) =(A(4)*A(8)-A(7)*A(5));

B(2,1) =-(A(2)*A(9)-A(8)*A(3));
B(2,2) =(A(1)*A(9)-A(7)*A(3));
B(2,3) =-(A(1)*A(8)-A(7)*A(2));

B(3,1) =(A(2)*A(6)-A(5)*A(3));
B(3,2) =-(A(1)*A(6)-A(4)*A(3));
B(3,3) =(A(1)*A(5)-A(4)*A(2));
% det(A)*inv(A)-B
end