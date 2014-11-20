% Compute the restoring force vectors
% by which the individual gcells act on the nodes.
% Return an array of elevec objects so they may be assembled.
%    Call as
% evs = restoring_force(feb, geom, u1, u), or
% [feb,evs] = restoring_force(feb, geom, u1, u)
%     geom     - geometry field
%     u1       - displacement field at the end of time step t_n+1
%     u        - displacement field at the end of time step t_n
%
% Note: This method *UPDATES* the state of the feblock object.  In
%       particular, the material state gets updated.  If this gets
%       called for converged u1, the feblock must be assigned to itself
%       on return; otherwise the first return value must be ignored!
%
function [self,evs] = restoring_force(self, geom, u1, u, alpha)
    evs= [];
end