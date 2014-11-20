function [lambda,v,converged]=gepbinvpwr2(K,M,v,tol,maxiter)
% Block inverse power method for the generalized eigenvalue problem.
%
% function [lambda,v,converged]=gepbinvpwr2(K,M,v,tol,maxiter)
%
% Block inverse power method for k smallest eigenvalues of the generalized
% eigenvalue problem
%           K*v= lambda*M*v
% 
% K,M = square matrix,
% v= initial guess of the eigenvectors (for instance random),
% tol= relative tolerance on the eigenvalue,
% maxiter= maximum number of allowed iterations
% Returns
% lambda = computed eigenvalue,
% v= computed eigenvector,
% converged= Boolean flag, converged or not?
% 
% (C) 2008, Petr Krysl
    [m,n] = size(K);
    if (m ~= n)
        error('Error: matrix K is not square');
    end
    [m,n] = size(M);
    if (m ~= n)
        error('Error: matrix M is not square');
    end
    nvecs=size(v,2);% How many eigenvalues?
    plambda=Inf+zeros(nvecs,1);% previous eigenvalue
    lambda =plambda;
    [L,U,p] =lu(K,'vector');
    converged = false;% not yet
    for iter=1:maxiter
        u=U\(L\(M*v(p,:))); % update vector
        for j=1:nvecs
            lambda(j)=(v(:,j)'*K*v(:,j))/(v(:,j)'*M*v(:,j));% Rayleigh quotient
        end
        [v,r]=qr(u,0);% economy factorization
norm(lambda-plambda)/norm(lambda)
        if (norm(lambda-plambda)/norm(lambda)<tol)
            converged = true; break;
        end
        plambda=lambda;
    end
end
