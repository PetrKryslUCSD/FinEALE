function midstep_F=compute_midstep_def_grad(F,prev_F)
% Compute the mid-step deformation gradient.
%
%    function midstep_F=compute_midstep_def_grad(F,prev_F)
%
  midstep_F= (F+prev_F)/2;
end
