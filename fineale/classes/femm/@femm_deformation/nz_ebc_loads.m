function F = nz_ebc_loads(self, assembler, geom, u)
 % Compute load vector for nonzero essential boundary conditions.
%
% function F = nz_ebc_loads(self, assembler, geom, u)
% 
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%
 error('Base class does not implement this behaviour!');
end