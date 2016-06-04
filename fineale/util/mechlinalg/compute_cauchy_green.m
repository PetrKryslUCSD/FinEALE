function C=compute_cauchy_green(F)
% Compute the Cauchy-Green deformation tensor.
%
% function C=compute_cauchy_green(F)
 C= F'*F;
end
