function out = cumspace(start,finish,spacings)
% Generate cumulative space.
% 
% function out = cumspace(start,finish,N)
%
% Generate a cumulative sequence of numbers between start and finish.
% This sequence corresponds to separation of adjacent numbers given by the
% vector 'spacings', stretched out or shrunk to fit between start and
% finish.

x=[0,cumsum(spacings)];
x =x/max(x);
out = start*(1-x)+finish*x;
end

