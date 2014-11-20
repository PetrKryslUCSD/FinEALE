% femm
%
% Classes for finite element model machines (femms).
% 
% Finite element model machine is an object that groups together services and 
% operations for a particular finite element model.
% For a given finite element model  one could create several finite element model 
% machines  (FEMMs).  For instance, for a single region domain one would have one 
% FEMM for the interior, and one FEMM  for each piece of the boundary with 
% prescribed nonzero conditions (nonzero flux, fixed nonzero temperature 
% or temperature rate).
%
% <a href="matlab: doc femm_base">femm_base</a> - Base class for finite element models
% <a href="matlab: doc femm_deformation_linear">femm_deformation_linear</a> - Linear mechanical deformation model
% <a href="matlab: doc femm_deformation_linear_sri">femm_deformation_linear_sri</a> - Linear mechanical deformation model with  selective reduced integration
% <a href="matlab: doc femm_heat_diffusion">femm_heat_diffusion</a> - Heat conduction model