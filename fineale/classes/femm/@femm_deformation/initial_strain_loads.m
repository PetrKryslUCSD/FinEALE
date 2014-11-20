% Compute the initial-strain load vectors of the individual gcells.
%
% function evs = initial_strain_loads(self, geom, u)
%
% Return an array of the element vectors so they may be assembled.
%    Call as
% ems = initial_strain_loads(femm, geom, u)
%     geom=geometry field
%     u=displacement field
%
function evs = initial_strain_loads(self, geom, u)
 error('Base class does not implement this behaviour!');
end