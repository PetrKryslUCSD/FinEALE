function [K,F]=apply_penalty_mpc(nfreedofs,dofnums,dofmultipliers,penfact)
% Apply multi point constraints (MPC) through the penalty technique.
%
% function [K,F]=apply_penalty_mpc(nfreedofs,dofnums,dofmultipliers,rhsvalue,penfact)
%
% Output:
% K= stiffness matrix
% F= global load vector 
%
% Input:
% nfreedofs= number of free degrees of freedom,
% dofnums= array of degree of freedom numbers,
% dofmultipliers= degree of freedom multipliers,
% penfact= penalty factor

    sv=zeros(nfreedofs,1);
    sv(dofnums)=reshape(dofmultipliers,length(dofnums),[]);
    K=penfact*sv*sv';
    F=penfact*sv;
end
