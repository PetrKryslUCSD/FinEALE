%% Single truss bar structure: stiffness and thermal load
%
% The structure is a single truss bar between points [10,-5,20] in and
% [30,25,-15] in. The bar is of uniform cross-section of 2 in. squared. The
% Young's modulus is assumed at value of 30e6 psi.  The structure is
% exposed to an increase of temperature of 100°F. The coefficient of
% thermal expansion is 5.5e-6 1/degree Fahrenheit. The goal is to calculate
% the stiffness matrix of the structure and the thermal load.

%% 
% The solution calculated in longhand is given by Rao, The Finite Element
% Method in Engineering, fifth edition, 2011, pages 315-317.

%% Solution

function pub_single_bar
pu= physical_units_struct;
%% 
% Define the finite element node set.
fens=fenode_set(struct('xyz',[...
    10,-5,20; ...% node 1
    30,25,-15 ...% node 2
    ]*pu.IN));

%% 
% Define the single finite element set. Note that we are setting the
% cross-sectional area as the property of the element set.
fes= fe_set_L2(struct('conn',[1,2], 'other_dimension', 2*pu.IN^2));


%% 
% Define the material property object and the material object. Note well
% the value of the Poisson's ratio  of zero.   Strictly speaking this is
% not necessary, the material stiffness calculation would correctly deduce
% the material  stiffness constant as the Young's modulus, but since it's
% value is not given in the reference we may take it as well as zero.
prop = property_deformation_linear_iso (struct (...
    'E',30e6*pu.PSI,...% Young's modulus
    'nu',0.,...
    'alpha', 5.5e-6*pu.IN/pu.IN/pu.F...% coefficient of thermal expansion
    ));
%% 
% Note that the material  object needs to be uniaxial in order to provide
% the correct calculation of the material stiffness matrix  (one by one),
% and the thermal stress (a single component, along the bar)
    mater = material_deformation_linear_uniax (struct('property',prop));
%% 
% Define the finite element model machine of the small-deformation linear
% elastic kind. One-point Gauss quadrature is in fact sufficient for the
% stiffness matrix and the thermal-strain load vector.
%% 
% Note well that we are making sure to define the material orientation
% matrix.   The truss bar has a special direction,  its axis.  The material
% orientation matrix  is a single vector  in the direction of the axis of
% the bar, and  is used to calculate the component of the displacement in
% the direction of the bar  axis and the strain-displacement matrix that
% produces axial strain from the  global 3-D displacements of the nodes. In
% this case the material orientation matrix  is provided by the finite
% element model machine itself, as requested by specifying the material
% orientation matrix as 'geniso'.
    femm = femm_deformation_linear(struct ('material',mater,...
        'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',1,'order',1)),...
        'Rm','geniso'));
%% 
% Define the geometry and displacement fields.
% Geometry is constructed from the finite element node set.
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
%% 
% The displacement field is created by cloning  the geometry field since
% they have the same number of parameters per node (3).
u   = 0*clone(geom,'u');
%% 
% The structure is free-floating: There are no EBC's applied on the
% structure. Therefore we proceed immediately to number equations:
u   = numberdofs (u);
%% 
% Finally we produce the first result: we assemble the element  (structure)  stiffness matrix:

K = stiffness(femm, sysmat_assembler_sparse, geom, u);

%% 
%  We print out the resulting stiffness matrix in the desired units:
format short e
full(K/(pu.LBF/pu.IN))


%% 
% And we may compare with the reference (hand-calculated)  stiffness matrix
% (in units of pounds per inch)..
ref_K=   1.0e+05 *[
    1.8916    2.8373   -3.3102   -1.8916   -2.8373    3.3102
    2.8373    4.2560   -4.9653   -2.8373   -4.2560    4.9653
   -3.3102   -4.9653    5.7929    3.3102    4.9653   -5.7929
   -1.8916   -2.8373    3.3102    1.8916    2.8373   -3.3102
   -2.8373   -4.2560    4.9653    2.8373    4.2560   -4.9653
    3.3102    4.9653   -5.7929   -3.3102   -4.9653    5.7929];
%% 
% Our calculation  is compared  to the reference stiffness matrix using matrix norm:
ekn=norm(full(K/(pu.LBF/pu.IN))- ref_K)/norm(ref_K)

%% 
% The accuracy of the calculation is actually better than suggested by the
% relative norm of the difference.  Only five digits are given for the reference
% matrix  and hence the difference  is relatively large.


%% 
% We calculate the thermal load vector. For this we need to define the nodal
% field that represents the temperature increment  across the structure.
dT = nodal_field(struct ('name',['dT'], 'dim', 1, 'nfens',fens.count));
%% 
% Note that it is necessary to use addition of the (originally zero) values
% in the  temperature field; we may also use the expression  |dT.values = dT.values + 100*pu.F|.
dT.values(:) = 100*pu.F;

%% 
% The thermal-strain  load vector for the element is assembled next....
F = thermal_strain_loads(femm, sysvec_assembler, geom, u, dT);
%% 
% ...  and the result is printed out  in pounds:
F/(pu.LBF)

%% 
% Reference (hand calculated) thermal load vector in pounds
ref_F=  1.0e+04 *[
  -1.313449091077187
  -1.970173636615779
   2.298535909385076
   1.313449091077187
   1.970173636615779
  -2.298535909385076];
%% 
% may be then compared with the calculated load vector.
efn=norm(F/(pu.LBF)-ref_F)/norm(ref_F)
%% 
% The norm  of the difference is close to machine epsilon, as in this case
% the reference vector is known to all significant digits.

%% Discussion
% 
%% 
% The toolbox provides answers identical to those calculated using
% elementary structural-analysis formulas for truss elements.