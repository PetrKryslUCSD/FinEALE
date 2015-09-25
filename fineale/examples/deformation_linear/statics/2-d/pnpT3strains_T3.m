% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Strain patterns calculated for a single triangle.
function pnpT3strains
E= 1.0; nu=  0.3;
%% 
xall= [-1, -1/2; 3, 2; 1, 2];%Coordinates of the nodes

conn= [1,2,3];% The definition of the element, listing its nodes

fens=fenode_set;
fens.xyz= xall;
fes=fe_set_T3(struct('conn',conn));
% Material
        prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
        mater = material_deformation_linear_biax (struct('property',prop, ...
            'reduction','stress'));
        % Finite element block
        femm = femm_deformation_linear(struct ('material',mater, 'fes',fes,...
            'integration_rule',tri_rule (struct('npts',1))));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 2, 'fens',fens));
        % Define the displacement field
        u   = clone(geom,'u');
        u   = u*0; % zero out
        
        % Number equations
        u   = numberdofs (u);
        % Assemble the system matrix
        K = stiffness(femm, sysmat_assembler_sparse,    geom, u);
        % Load