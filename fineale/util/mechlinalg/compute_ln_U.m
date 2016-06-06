function lnU=compute_ln_U (U)
%   Compute the logarithmic strain tensor.
%
%  function lnU=compute_ln_U (U)
% 
%   Given the U tensor, compute its logarithm.
%   An eigendecomposition is needed to perform
%   this task.

[evecs,evals]=eig(U);
lnevals = log(diag(evals));
lnU=evecs*diag(lnevals)*evecs';
end
