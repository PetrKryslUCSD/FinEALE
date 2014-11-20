function Rfield= compose_rotation_fields(RfieldA, RfieldB)
% Compose two rotation fields.
%
% function Rfield= compose_rotation_fields(RfieldA, RfieldB)
%
Rfield=RfieldA;
RsA = RfieldA.values;
RsB = RfieldB.values;
for i=1:size(RsA,1)
    RA= reshape(RsA(i,:),3,3);
    RB= reshape(RsB(i,:),3,3);
    RsA(i,:)= reshape(RB*RA,1,9);
end
Rfield.values= RsA;
end

