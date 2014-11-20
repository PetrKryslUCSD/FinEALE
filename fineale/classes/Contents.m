% Each folder corresponds to a class tree. The folders whose names
% start with the @ each hold the methods defined for the class. For
% instance, @fe_set_H20 holds the methods of the fe_set_H20 class.    
%
% Note:
% fineale objects are derived as Matlab classes. 
% Refer to the <a href="matlab:doc 'object-oriented-design-with-matlab'">Matlab documentation of  the object-oriented features.</a>
% 
%
% Class trees:
%
% <a href="matlab: doc 'fineale/classes/data_colormap/Contents.m'">data_colormap</a> - Mapping of values to colors. 
% <a href="matlab: doc 'fineale/classes/fe_set/Contents.m'">fe_set</a> - Group of classes to represent finite elements.
% <a href="matlab: doc 'fineale/classes/femm/Contents.m'">femm</a> - Group of classes to represent finite element model machines.
% <a href="matlab: doc 'fineale/classes/fenode_set/Contents.m'">fenode_set</a> - Class for finite element nodes.
% <a href="matlab: doc 'fineale/classes/fenode_to_fe_map/Contents.m'">fenode_to_fe_map</a> - Class for  mapping of nodes to finite elements.
% <a href="matlab: doc 'fineale/classes/force_intensity/Contents.m'">force_intensity</a> - Class for representation of distributed loads.
% <a href="matlab: doc 'fineale/classes/graphic_viewer/Contents.m'">graphic_viewer</a> - Graphic viewer class.
% <a href="matlab: doc 'fineale/classes/integration_rule/Contents.m'">integration_rule</a> - Group of classes for numerical integration rules
% <a href="matlab: doc 'fineale/classes/material/Contents.m'">material</a> - Group of classes for operations on materials.
% <a href="matlab: doc 'fineale/classes/nodal_field/Contents.m'">nodal_field</a> - Class for nodal fields  as representations of variables on the mesh.
% <a href="matlab: doc 'fineale/classes/property/Contents.m'">property</a> - Group of classes for material properties. 
% <a href="matlab: doc 'fineale/classes/sysmat_assembler/Contents.m'">sysmat_assembler</a> - Group of classes for assembling of  system matrices.
% <a href="matlab: doc 'fineale/classes/sysvec_assembler/Contents.m'">sysvec_assembler</a> - Class for assembling of system vectors.
%
%
% Notes: 
% 1. One may find out which attributes are defined for a class by
%    typing the constructor on the command line.
%    For instance, executing `fe_set_H20' will print a list of attributes 
%    that can be obtained from the object and those that can be set.
% 2. One may inquire which methods are defined for a class
%    using the `methods' Matlab command. For instance, `methods nodal_field'.
% 3. Constructors always take either zero or one argument.  For no
%    arguments, the default object is created.  For one argument, the argument
%    may be either an object of the same type, in which case a copy is
%    created; or, it may be a struct with one or more fields
%    representing the parameters to be passed to the constructor. 
%      For instance, these are equivalent calls to the
%    finite element node constructor: 
%       f=fenode_set(struct('xyz', [0, -.2])) 
%       Parameters.xyz=[0, -.2]; f=fenode_set(Parameters)
% 4. Attributes may be retrieved from or changed in an object using the
%    "dot" access methods. For instance, to retrieve the coordinates from 
%    the finite element node set  use
%         Parameters.id=5; Parameters.xyz=[0, -.2]; f=fenode_set(Parameters)
%         f =
%           fenode_set
%
%           Properties:
%             xyz: [0 -0.2000]
%           Methods
%         f.xyz
%         ans =
%                  0   -0.2000
%    To modify the coordinates in the finite element node set use
%         f.xyz= [33, 13]
%         f =
%           fenode_set
%
%           Properties:
%             xyz: [33 13]
%           Methods
% 5. The object for which a method gets called is always called "self" in 
%    the method signature.
