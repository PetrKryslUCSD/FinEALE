function Volumes =t4util_volume (t,v)
% Compute the volumes of tetrahedra.
%
% function Volumes =t4util_volume (t,v)
% 
    Volumes = 0*double(t(:,1));
    for iS1 =1:size(t,1)
         Volumes (iS1) =tetvol(v(t(iS1,:),:));
    end
end
function vol = tetvol(X)
        vol = det([X(2,:)-X(1,:);X(3,:)-X(1,:);X(4,:)-X(1,:)])/6;
    end