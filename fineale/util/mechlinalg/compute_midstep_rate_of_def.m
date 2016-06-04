function [midstep_D,midstep_W]=compute_midstep_rate_of_def (rel_F, dt)
% Compute the mid-step rate of deformation sensor.
%
% function [midstep_D,midstep_W]=compute_midstep_rate_of_def (rel_F, dt)
s= (eye(3)+ rel_F)/2;
q= rel_F-eye(3);
r=s\q;

midstep_D = 0.5/dt * (r + r');
midstep_W = 0.5/dt * (r - r');
end