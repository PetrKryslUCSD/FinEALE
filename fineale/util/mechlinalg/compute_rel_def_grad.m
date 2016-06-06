function rel_F=compute_rel_def_grad(F,prev_F)
% Computes the relative deformation gradient.
%
% function rel_F=compute_rel_def_grad(F,prev_F)
%
% F= Current configuration deformation gradient  (d x_n+1/d X),
% prev_F= Previous configuration deformation gradient (d x_n/d X)
%
% rel_F= relative deformation gradient,  (d x_n+1/d x_n)

det_prev_F = det(prev_F);
if (det_prev_F <= 0)
    error('Failed to invert def grad');
end

rel_F= F/prev_F;
end
