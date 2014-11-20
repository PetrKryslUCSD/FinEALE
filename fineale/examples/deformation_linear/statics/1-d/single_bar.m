
function single_bar
% Single truss bar calculation.
pu= physical_units_struct;
%% 
% Define the finite element nodes.
fens=fenode_set(struct('xyz',[...
    10,-5,20; ...% node 1
    30,25,-15 ...% node 2
    ]*pu.IN));

%% 
% Define the single finite element.
fes= fe_set_L2(struct('conn',[1,2], 'other_dimension', 2*pu.IN^2));


%% 
% Define the material property object and the material object.
prop = property_deformation_linear_iso (struct ('E',30e6*pu.PSI,'nu',0));
    mater = material_deformation_linear_uniax (struct('property',prop));
%% 
% Define the finite element model machine of the small-deformation linear
% elastic kind. One-point Gauss quadrature would be sufficient for the
% stiffness matrix, but the load needs the two point quadrature.
    femm = femm_deformation_linear(struct ('material',mater,...
        'fes',fes,...
        'integration_rule',gauss_rule(struct('dim',1,'order',2)),...
        'Rm','geniso'));
%% 
% Define the geometry and displacement fields.
% Geometry is constructed from the finite element node set.
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
%% 
% The displacement field is created by cloning  the geometry field since
% they have the same number of parameters per node (3).
u   = 0*clone(geom,'u');
% There are no EBC's applied on the structure.

% Number equations
u   = numberdofs (u);
% Assemble the system matrix
K = stiffness(femm, sysmat_assembler_sparse, geom, u);
full(K/(pu.LBF/pu.IN))


%% 
% Reference (hand-calculated)  stiffness matrix in pounds per inch
ref_K=   1.0e+05 *[
    1.8916    2.8373   -3.3102   -1.8916   -2.8373    3.3102
    2.8373    4.2560   -4.9653   -2.8373   -4.2560    4.9653
   -3.3102   -4.9653    5.7929    3.3102    4.9653   -5.7929
   -1.8916   -2.8373    3.3102    1.8916    2.8373   -3.3102
   -2.8373   -4.2560    4.9653    2.8373    4.2560   -4.9653
    3.3102    4.9653   -5.7929   -3.3102   -4.9653    5.7929];
norm(full(K/(pu.LBF/pu.IN))- ref_K)/norm(ref_K)
% Load
