function K = stiffness (self, assembler, geom, u)
% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, u)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u=displacement field
 error('Base class does not implement this behaviour!');
end