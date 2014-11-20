% Perform isoparametric interpolation in 1D
% 
% function v=isoparamipol1d(gcells,pc,fldv)
% 
% gcell= geometric cell
% pc= array of parametric coordinates
% fldv= field values at nodes of the cell
% 
function v=isoparamipol1d(gcell,pc,fldv)
v=zeros(1,length(pc));
for j=1:length(pc)
    N=bfun(gcell,pc(j));
    v(j)=N'*fldv;
end 