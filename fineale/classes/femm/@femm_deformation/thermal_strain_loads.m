function F = thermal_strain_loads(self, assembler, geom, u, dT)
% Compute the thermal-strain load vectors of the finite element set.
%
% function F = thermal_strain_loads(self, assembler, geom, u, dT)
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%     dT=temperature difference field (current temperature minus the
%         reference temperature at which the solid experiences no strains)
%
 error('Base class does not implement this behaviour!');
end