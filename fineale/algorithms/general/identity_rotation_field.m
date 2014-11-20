function Rfield0= identity_rotation_field(name, nfens)
% Create an identity rotation field.
%
% function Rfield0= identity_rotation_field(name, nfens)
%
% Example:
% Rfield0= identity_rotation_field('Rfield', geom0.nfens);

Rfield0 = nodal_field(struct ('name',[name], 'dim', 9, 'data', ones(nfens,1)*reshape(eye(3),1,9)));
end

