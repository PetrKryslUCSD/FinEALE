function [K,F]=apply_penalty_mpc(nfreedofs,dofnums,umultipliers,rhsvalue,penfact)
% Apply multi point constraints (MPC) through the penalty technique.
%
% function [K,F]=apply_penalty_mpc(nfreedofs,dofnums,umultipliers,rhsvalue,penfact)
%
% K= stiffness matrix
% F= global load vector 
% uebc = field which describes the constraints,
% u= field which does not have the constraints applied, and serves as the source of equation numbers,
% penfact= penalty multiplier
    sv=zeros(nfreedofs,1);
    sv(dofnums)=reshape(umultipliers,length(dofnums),[]);
    K=penfact*sv*sv';
    F=penfact*sv;
end
