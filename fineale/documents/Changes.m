% 04/11/2015
% - Abaqus export added from the research files.
% - Added example with  2-D surface embedded in 3-D.

% 03/25/2015
% - Sparse system matrix assembly reconfigured to use the trick with an additional unknown.
% The assembly is speeded up by about 30%.
% - Nodal field: added a method for extracting the degree of freedom numbers in vector form
% in one shot instead of using two methods.  Faster.
% - drawmesh(): Create rotation field for use by beam elements.

% 01/01/2015
% - mirror_mesh(): was  not using Point -- corrected.

% 12/18/2014
% - Rm use in heat  conduction: Rm is now used to produce the "reduced"
% Jacobian matrix.

% 12/13/2014
% - merge_nodes()  was corrected to merge multiple  nodes within the same
% location.  Previously only two nodes  were merged.
% - fuse_nodes() now chooses from the same-location set the lowest number of node.

% 12/09/2014
% - Added passing of function handles to  set  values of  essential boundary conditions.
% Affected algorithms: deformation_linear_statics, heat_diffusion_steady_state

% 12/08/2014
% - fe_set: A check in the setter function to replace  connectivity is now
% used to compare the number of finite element nodes with  the value of the
% field

% 12/07/2014
% - mesh_boundary(): An additional algorithm has been added to compensate for missing
% functionality in Julia. The operation of the algorithm  (the results)
% do not change.
% - The clear() function of the graphic_viewer class does not clear the figure
%     only the axes.

% 11/12/2014
% - deformation_plot_deformation().   Additional options.

% 11/08/2014
% - Added nonlinear benchmark problems.

% 09/13/2014
% - Added  nonlinear static deformation algorithm.
% - New method is used for changing the state of an FEMM.
%   Instead of restoring_forces(), the update() method is now used.

% 06/11/2012
% - Added  implicit time stepping for linear dynamics. Updated documentation.

% 06/10/2012
% - Added  explicit time stepping for linear dynamics.

% 05/27/2012
% – Added factorization to transient heat conduction and modal analysis.
% – Updated documentation of algorithms.

% 05/26/2012
% - Added mesh generation examples from FAESOR.
% - Cleaned up notation in the classes.
% - Redesigned renumbering of the  equations.   Note that the nodal_field method
% numbereqns() had to change incompatibly with the respect to FAESOR.

% 05/17/2012
% - deformation_linear_statics() description of the model_data input needs to be provided.
% - for constant material stiffness the recomputation of the tangent moduli should be avoided.

% 05/12/2012
% – The heat source and flux load are now handled using the distrib_loads method.
 
% 05/10/2012
% – Implemented algorithms for heat conduction. The example files
% may now be rewritten using algorithms which should lead to easier maintenance.

% 05/08/2012
% – First workable implementation  of the classes based on  
% the new Matlab classdef mechanism.

 