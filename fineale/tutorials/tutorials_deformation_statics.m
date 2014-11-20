% Static stress analysis
%
% <a href="matlab:showdemo 'pub_single_bar'">Single truss bar structure: stiffness and thermal load.</a> 
%     Features:
%     - Static stress analysis of 3-D truss member. 
%     - Stiffness matrix and thermal-strain load vector compared with
%       analytical results.
%     - Use of material orientation matrix to formulate uniaxial model of deformation.
%
% <a href="matlab:showdemo 'pub_plane_w_hole'">Infinite plane with a circular hole.</a> 
%     Features:
%     - Static stress analysis of equivalent 2-D and 3-D models. 
%     - Stress concentration factor compared with analytical results.
%     - Application of exact stress components as boundary conditions.
%
% <a href="matlab:showdemo 'pub_LE10NAFEMS_T10'">Elliptical plate with elliptical hole; T10 tetrahedra.</a> 
%     Features:
%     - Static stress analysis of true 3-D geometry  with curved surfaces.
%     - Meshed with quadratic T10 tetrahedral elements.
%     - Extraction of stress by extrapolation from quadrature points.
%     - Use of Superconvergent Patch Recovery (SPR).
%
% <a href="matlab:showdemo 'pub_LE10NAFEMS_H64'">Elliptical plate with elliptical hole; H64 elements with nodal quadrature.</a> 
%     Features:
%     - Static stress analysis of true 3-D geometry  with curved surfaces.
%     - Meshed with cubic H64 hexahedral elements.
%     - Use of nodal quadrature (Newton-Cotes rule).
%
% <a href="matlab:showdemo 'pub_Floyd'">Floyd's pressure vessel.</a> 
%     Features:
%     - Static stress analysis of axially symmetric geometry.
%     - Automatically meshed with quadratic triangular elements.
%     - Calculation of principal stresses.
%     - Interpolation of nodal field at other points within the mesh.
%
% <a href="matlab:showdemo 'pub_thick_pipe'">Thick pipe with internal pressure.</a> 
%     Features:
%     - Static stress analysis of pressurized pipe geometry with analytical 
%       solution.
%     - Meshed with tetrahedral and hexahedral elements.
%     - Use of reduced integration, uniform and selective, to treat nearly 
%       incompressible materials.
%     - Access to integration-point data.
%
% <a href="matlab:showdemo 'pub_thick_pipe_ps'">Thick pipe with internal pressure: plane strain model.</a> 
%     Features:
%     - Static stress analysis of pressurized pipe geometry with analytical 
%       solution.
%     - Uses the plane-strain model dimension reduction.
%     - Use of reduced integration, uniform and selective, to treat nearly 
%       incompressible materials.
%     - Access to integration-point data using output orientation matrix.
%
% <a href="matlab:showdemo 'pub_thick_pipe_axi'">Thick pipe with internal pressure: axially symmetric model.</a> 
%     Features:
%     - Static stress analysis of pressurized pipe geometry with analytical 
%       solution.
%     - Uses the axially symmetric model dimension reduction.
%     - Use of reduced integration, uniform and selective, to treat nearly 
%       incompressible materials.
%     - Access to integration-point data.
%
% <a href="matlab:showdemo 'pub_R0031NAFEMS'">Laminated Strip Under Three-Point Bending.</a> 
%     Features:
%     - Static stress analysis of laminated plate.
%     - Meshing of layered plates.
%     - One region per set of layers with the same orientation.
%     - Calculation  of stress components at various points in the layup.
%
% <a href="matlab:showdemo 'pub_R0031NAFEMS_1_region'">Laminated Strip Under Three-Point Bending: a single-region formulation.</a> 
%     Features:
%     - Static stress analysis of laminated plate.
%     - Meshing of layered plates.
%     - Calculation  of stress components at various points in the layup 
%       based on a single region with material orientation function.
%
% <a href="matlab:showdemo 'pub_Twisted_beam'">Twisted beam  benchmark.</a> 
%     Features:
%     - Static stress analysis of a twisted beam (Macneal, Harder benchmark).
%     - Use of Selective Reduced Integration (SRI) to improve bending 
%       response in hexahedra.
%     - Convergence of displacement, running simulation with progressive
%       refinement.
%
% <a href="matlab:showdemo 'pub_rodwtrans'">Stress  calculation in a tension rod with shoulder fillet.</a> 
%     Features:
%     - Static axially-symmetric stress analysis of rod under tensile load.
%     - Variety of boundary conditions considered.
%     - A variation of stress highlighted with elevation plots.
%