function lnU=compute_ln_U (U)
%   Compute the logarithmic strain tensor.
%   Given the U tensor, compute its logarithm.
%   An eigendecomposition is needed to perform
%   this task.
%
%  function lnU=compute_ln_U (U)
% 

[evecs,evals]=eig(U);
lnevals = log(diag(evals));
lnU=evecs'*diag(lnevals)*evecs;
end
