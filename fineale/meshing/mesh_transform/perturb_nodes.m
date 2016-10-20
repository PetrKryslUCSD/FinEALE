function [fens] = perturb_nodes(fens, PerturbationAmplitudeVector, free)
% Perturb locations of nodes by a random amount.
%
% [fens] = perturb_nodes(fens, PerturbationAmplitudeVector, free)
%
% fens= Finite element nodes of the mesh, 
% PerturbationAmplitudeVector= m-D vector where m is the number of spatial
%      coordinates
% free  = optional: array of ones  and zeros.  One row per node.
%      In the columns: 1 if the note is free to move in this direction, 0
%      if the node  is not allowed to move in this direction. 
% 
% 
if (~exist('free', 'var'))
    free= ones(size(fens.xyz));
end
r=(rand(size(free))- 0.5)*2;
for   j=1:size(fens.xyz,2)
    fens.xyz(:,j)=fens.xyz(:,j)+(r(:,j).*free(:,j))*PerturbationAmplitudeVector(j);
end
end
