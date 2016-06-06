function [midstep_D,midstep_W]=compute_midstep_rate_of_def (rel_F, dt)
% Compute the mid-step rate of deformation sensor.
%
% function [midstep_D,midstep_W]=compute_midstep_rate_of_def (rel_F, dt)
%
% rel_F= relative deformation gradient  (d x_n+1/d x_n)
%
% rel_F= relative deformation gradient,  (d x_n+1/d x_n)
% dt= time step

rel_F_mid= (eye(3)+ rel_F)/2;
du_dxn= rel_F-eye(3);
h=du_dxn/rel_F_mid;% Relative incremental displacement gradient

midstep_D = 0.5/dt * (h + h');
midstep_W = 0.5/dt * (h - h');
end