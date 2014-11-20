function [K,F]=apply_penalty_ebc(K,F,uebc,u,penfact)
% Apply penalty essential boundary conditions.
%
% function [K,F]=apply_penalty_ebc(K,F,uebc,u,penfact)
%
% K= stiffness matrix
% F= global load vector 
% uebc = field which describes the constraints,
% u= field which does not have the constraints applied, and serves as the source of equation numbers,
% penfact= penalty multiplier
    ip=gather(uebc,(1:get(uebc,'nfens')),'is_fixed','noreshape');
    ix=find (ip~=0);
    pv=gather_fixed_value(uebc,(1:get(uebc,'nfens')));
    dofnums=gather_dofnums(u,(1:get(u,'nfens')));
    dofnums=dofnums(ix);
    pv=pv(ix);
    penalty=penfact*max(diag(K));
    for j=1:length(dofnums)
        gj=dofnums(j);
        K(gj,gj) = K(gj,gj) + penalty;
    end
    for j=1:length(dofnums)
        gj=dofnums(j);
        F(gj) = F(gj) + penalty*pv(j);
    end
end
