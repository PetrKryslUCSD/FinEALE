% Compute the mass matrix by computing and assembling the
% matrices of the individual FEs.
%
% function M = mass (self, assembler, geom, u)
%
% Return a mass matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u=displacement field
function M = mass (self, assembler, geom, u)
M = mass@femm_deformation_linear(self, assembler, geom, u);
M(find(diag(M)<0))=0;% This is a horrible kludge  intended to fix problems with negative-volume  elements
end


