function elerrs = flux_L2_error (self, geom, u, dT, nodal_stress)
% Compute the L2 stress flux error of the individual gcells.
%
% function elerrs = flux_L2_error (self, geom, u, dT, nodal_stress)
%
% Input arguments
%     self     - finite element block
%     geom     - reference geometry field
%     u        - displacement field
%     dT       - temperature difference field, may be supplied as empty
%     nodal_stress   -  field of the nodal values of the stress. This would be computed as
%       nodal_stress = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', 1:4);
%       for a plane strain problem (as an example).  
% 
% Output argument
% Return an array of the element L2 errors of the flux.
%     elerrs= Array, number of columns equals the number of geometric
%     cells; for each finite element e the array component elerrs(e) holds the 
%     L2 norm of the error of the flux over the element 
%      
  error('Base class does not implement this behaviour!');
end