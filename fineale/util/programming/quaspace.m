function out = quaspace(start,finish,N)
% Generate quadratic space.
% 
% function out = quaspace(start,finish,N)
%
% Generate a quadratic sequence of numbers between start and finish.
% This sequence corresponds to separation of adjacent numbers that
% increases linearly from start to finish.
x=linspace(0,1,N);
x=cumsum(x);
x=cumsum(x);
x =x/max(x);
out = start*(1-x)+finish*x;
end

