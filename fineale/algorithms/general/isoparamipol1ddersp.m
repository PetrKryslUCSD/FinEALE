function v=isoparamipol1ddersp(gcell,pc,xv,fldv)
% Perform isoparametric interpolation in 1D, derivatives.
% 
% function v=isoparamipol1ddersp(gcell,pc,xv,fldv)
% 
% gcell= geometric cell
% pc= array of parametric coordinates
% fldv= field values at nodes of the cell
% 
v=zeros(1,length(pc));
for j=1:length(pc)
    Nder = bfundpar (gcell,pc(j));
    Ndersp = Nder/(xv'* Nder);
    v(j)=Ndersp'*fldv;
end 