
classdef nodal_field
    % Class to represent nodal fields.
    %
    % Class that represents nodal fields: geometry, displacement, incremental rotation,
    % temperature,... defined by values associated with nodes:
    %   They are defined as the finite
    %   element expansion of the unknowns, i.e. a field of the form
    %   u_h(x) = sum_i N_i(x)*u_i
    %   where N_i(x) are the basis functions, and u_i are the
    %   unknown nodal parameters (degrees of freedom).
    %
    
    properties
        name = '';% Name of the field
        values = [];% Array of degree of freedom parameters,  indexed by node
        dofnums = [];% Array of degree of freedom numbers, indexed by node
        is_fixed = [];% Array of Boolean flags, indexed by node
        fixed_values = [];% Array of fixed values, indexed by node
        nfreedofs = 0;% Total number of free degrees of freedom
    end
    
    methods
        
        
        function self = nodal_field (Parameters)
        % Constructor.
        %
        %  Parameters:
        % always required:
        %     name=name of the field
        % In addition, supply one of the following combinations:
        %     dim=dimension of the nodal parameters
        %         (1=scalar field,
        %          2=2D vector field,
        %          3=3D vector field, etc)
        %     nfens=number of nodal parameters in this field (equal to the number
        %     of nodes)
        %   or
        %     dim=dimension of the nodal parameters
        %         (1=scalar field,
        %          2=2D vector field,
        %          3=3D vector field, etc)
        %     fens=array of objects of the class fenode
        %   or
        %     data=data array of the nodal parameters
        %
        %   The second form is useful to initialize a geometry field (i.e. a field
        %   where the parameters are the coordinates of the nodes). The third form
        %   is useful when creating a field out of of a data array, rows
        %   corresponding to nodal parameters, number of columns is the dimension
        %   of the field.
        %
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            class_name='field';
            if (nargin == 0)
                Parameters=struct( [] );
            end
            self.name = '';
            if (isfield(Parameters,'name'))
                self.name  = Parameters.name;;
            end
            if (isfield(Parameters,'fens'))% we are creating the geometry field
                dim  = Parameters.dim;
                x = Parameters.fens.xyz;
                nfens = size(x,1);
                self.values = zeros(nfens,dim);
                self.values = x;
            elseif (isfield(Parameters,'nfens'))% general field based on the number of nodes
                dim  = Parameters.dim;
                nfens = Parameters.nfens;
                self.values = zeros(nfens,dim);
            elseif (isfield(Parameters,'data'))% general field based on supplied data
                self.values = Parameters.data;
            end
            self.dofnums = zeros(size(self.values));
            self.is_fixed = zeros(size(self.values));
            self.fixed_values = zeros(size(self.values));
            self.nfreedofs = 0;
        end
        
    end
    
    methods
        
        function  val=reshape (self,val)
        % Change the shape (reshape) an attribute of the field.
        %
        % function  val=reshape (self,val)
        %
        % First call gather_xxx(), and then reshape() (if needed).
        %
        % The shape returned is a vector (column matrix). 
        %
        % Note: this needs to be the non-conjugate  transpose.  We are
        % simply interested the same numbers, except reshaped into a
        % different array.
            val = reshape(val.',[],1);
        end
        
        function val = gather_values (self, conn)
        % Gather the field degree of freedom values.
        %
        % function val = gather_values (self, conn)
        %
        % Gather the field values corresponding to the finite element node
        % indices.
        %
        %       conn  - node numbers for which to gather stuff; if supplied as
        %                 empty, all nodes are assumed to be implied; and
        %
        %       The actual data in rectangular array form are returned.
        %       If desired,    the data may be reshaped into a column vector of
        %       (n*dim) rows, where n=length(conn) and dim=dimension
        %       of the nodal parameter of the field, with the components ordered
        %       consecutively, using the  reshape() method.
        %
        %   For instance, if
        %       f.dofnums=[ 1 2 3; 4 0 5; 0 0 6 ]
        %   and conn=[1 3]
        %       gather_dofnums(f,conn) yields [ 1 2 3; 0 0 6 ]
        %
            if (isempty(conn))
                nfens = size(self.values, 1);
                conn = (1:nfens);
            end
            val = self.values(conn,:);
        end
        
        function val = gather_is_fixed (self, conn)
        % Gather the field degree of freedom Boolean flags indicating whether or not it is fixed.
        %
        % function val = gather_is_fixed (self, conn)
        %
        % Gather the field values corresponding to the finite element node
        % indices.
        %
        %       conn  - node numbers for which to gather stuff; if supplied as
        %                 empty, all nodes are assumed to be implied; and
        %
        %       The actual data in rectangular array form are returned.
        %       If desired,    the data may be reshaped into a column vector of
        %       (n*dim) rows, where n=length(conn) and dim=dimension
        %       of the nodal parameter of the field, with the components ordered
        %       consecutively, using the  reshape() method.
        %
        %   For instance, if
        %       f.dofnums=[ 1 2 3; 4 0 5; 0 0 6 ]
        %   and conn=[1 3]
        %       gather_dofnums(f,conn) yields [ 1 2 3; 0 0 6 ]
        %
            if (isempty(conn))
                nfens = size(self.values, 1);
                conn = (1:nfens);
            end
            val = self.is_fixed(conn,:);
        end
        
        
        function val = gather_dofnums (self, conn)
        % Gather the field degree of freedom numbers.
        %
        % function val = gather_dofnums (self, conn)
        %
        % Gather the field degree of freedom numbers corresponding to the finite element node
        % indices. Note: only the free degrees of freedom are given nonzero 
        % numbers; the fixed degrees of freedom are given the number zero (0).
        %
        %       conn  - node numbers for which to gather stuff; if supplied as
        %                 empty, all nodes are assumed to be implied; and
        %
        %       The actual data in rectangular array form are returned.
        %       If desired,    the data may be reshaped into a column vector of
        %       (n*dim) rows, where n=length(conn) and dim=dimension
        %       of the nodal parameter of the field, with the components ordered
        %       consecutively, using the  reshape() method.
        %
        %   For instance, if
        %       f.dofnums=[ 1 2 3; 4 0 5; 0 0 6 ]
        %   and conn=[1 3]
        %       gather_dofnums(f,conn) yields [ 1 2 3; 0 0 6 ]
        %
        if (isempty(conn))% return all
            val = self.dofnums;
        else
            val = self.dofnums(conn,:);
        end
        end
        
        function val = gather_dofnums_vec (self, conn)
        % Gather the field degree of freedom numbers as a vector.
        %
        % function val = gather_dofnums_vec (self, conn)
        %
        % Gather the field degree of freedom numbers corresponding to the finite element node
        % indices. Note: only the free degrees of freedom are given nonzero 
        % numbers; the fixed degrees of freedom are given the number zero (0).
        %
        %       conn  - node numbers for which to gather stuff; if supplied as
        %                 empty, all nodes are assumed to be implied; and
        %
        %       The data is returned in the form of a column vector.
        %
        %   For instance, if
        %       f.dofnums=[ 1 2 3; 4 0 5; 0 0 6 ]
        %   and conn=[1 3]
        %       gather_dofnums_vec(f,conn) yields [ 1 2 3 0 0 6 ]'
        %
        val = reshape( self.dofnums(conn,:).',[],1);
        end
        
        function val = gather_fixed_values (self, conn)
              % Gather the field fixed degree of freedom values.
        %
        % function val = gather_fixed_values (self, conn)
        %
        % Gather the field fixed values corresponding to the finite element node
        % indices.
        %
        %       conn  - node numbers for which to gather stuff; if supplied as
        %                 empty, all nodes are assumed to be implied; and
        %
        %       The actual data in rectangular array form are returned.
        %       If desired,    the data may be reshaped into a column vector of
        %       (n*dim) rows, where n=length(conn) and dim=dimension
        %       of the nodal parameter of the field, with the components ordered
        %       consecutively, using the  reshape() method.
        %
        %   For instance, if
        %       f.dofnums=[ 1 2 3; 4 0 5; 0 0 6 ]
        %   and conn=[1 3]
        %       gather_dofnums(f,conn) yields [ 1 2 3; 0 0 6 ]
        %
            if (isempty(conn))
                nfens = size(self.values, 1);
                conn = (1:nfens);
            end
            val = self.fixed_values(conn,:);
        end
        
        function val = dim(self)
        % Dimension of the nodal parameters (i. e.  how many degrees of freedom per node).
            val = size(self.values, 2);
        end
        
        function val = nfens(self)
        % Number of nodes associated with the field.
            val = size(self.values, 1);
        end
        
    end
    
end


