function Rfield= update_rotation_field(Rfield, Incremental_rotation)
% Update rotation field by incremental rotation. 
% 
% function Rfield= update_rotation_field(Rfield, Incremental_rotation)
% 
    Rs = Rfield.values; 
    dA = Incremental_rotation.values; 
    for i=1:size(Rs,1)
        R= reshape(Rs(i,:),3,3);
        Rs(i,:)= reshape(rotmat(dA(i,:))*R,1,9);
    end 
    Rfield.values= Rs;
end

