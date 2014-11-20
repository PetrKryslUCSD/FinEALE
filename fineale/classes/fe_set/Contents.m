% fe_set
%
%   A fe_set is an object for a set of finite elements of the same type.  
%   The fe_set class is abstract (i.e. there is no point 
%   in making objects of this class), but all other fe_sets
%   (fe_set_H8 for hexahedra, ...) are descendants of this 
%   class. More precisely, they are descendents of the classes
%   fe_set_X_manifold, which endow the class fe_set with a
%   dimension (0 for a point, 1 for a curve, and so on).
% 
% <a href="matlab: doc fe_set">fe_set class documentation</a> 
%
% <a href="matlab: doc fe_set_0_manifold">fe_set_0_manifold class documentation</a> 
%
% <a href="matlab: doc fe_set_P1">fe_set_P1 class documentation</a> 
%
% <a href="matlab: doc fe_set_1_manifold">fe_set_1_manifold class documentation</a> 
%
% <a href="matlab: doc fe_set_L2">fe_set_L2 class documentation</a> 
% <a href="matlab: doc fe_set_L3">fe_set_L3 class documentation</a> 
% <a href="matlab: doc fe_set_L4">fe_set_L4 class documentation</a> 
%
%
% <a href="matlab: doc fe_set_2_manifold">fe_set_2_manifold class documentation</a> 
%
% <a href="matlab: doc fe_set_T3">fe_set_T3 class documentation</a> 
% <a href="matlab: doc fe_set_T6">fe_set_T6 class documentation</a> 
% <a href="matlab: doc fe_set_Q4">fe_set_Q4 class documentation</a> 
% <a href="matlab: doc fe_set_Q8">fe_set_Q8 class documentation</a> 
% <a href="matlab: doc fe_set_Q9">fe_set_Q9 class documentation</a> 
% <a href="matlab: doc fe_set_Q16">fe_set_Q16 class documentation</a> 
%
% <a href="matlab: doc fe_set_3_manifold">fe_set_3_manifold class documentation</a> 
% 
% <a href="matlab: doc fe_set_T4">fe_set_T4 class documentation</a> 
% <a href="matlab: doc fe_set_T10">fe_set_T10 class documentation</a> 
% <a href="matlab: doc fe_set_H8">fe_set_H8 class documentation</a> 
% <a href="matlab: doc fe_set_H20">fe_set_H20 class documentation</a> 
% <a href="matlab: doc fe_set_H27">fe_set_H27 class documentation</a> 
% <a href="matlab: doc fe_set_H64">fe_set_H64 class documentation</a> 
