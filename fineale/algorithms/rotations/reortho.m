% Re-orthogonalize an approximate rotation matrix
%
% function Rn=reortho(R)
%
function Rn=reortho(R)
Rn=R;
Rn(:,1)=Rn(:,1)/norm(Rn(:,1));
Rn(:,3)=skewmat(Rn(:,1))*Rn(:,2);
Rn(:,3)=Rn(:,3)/norm(Rn(:,3));
Rn(:,2)=skewmat(Rn(:,3))*Rn(:,1);    
Rn(:,2)=Rn(:,2)/norm(Rn(:,2));
