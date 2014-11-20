classdef sysmat_assembler_base < handle
% The class sysmat_assembler is for assembling of a sparse matrix from elementwise matrices.
% This is an abstract class, only its subclasses can be used.  This defines just an interface.
%
% Note that this class is derived from the handle class.   The properties of the 
% class may be changed from the outside –- the variables of this class are passed 
% by reference not by value.
%
% All descendents of this class  must implement the  following methods:
% function start_assembly(self, elem_mat_nrows, elem_mat_ncols, elem_mat_nmatrices, ndofs_row, ndofs_col)
% function assemble(self, mat, dofnums_row, dofnums_col)
% function assemble_symmetric(self, mat,  nums)
% function S= make_matrix (self)
%
	properties (Constant, GetAccess = public)
    % Fixed degrees of freedom numbers are given this value:
     % it indicates that this is not a valid free  degree of freedom number.
		fixed_dofnum=0;
	end
    
	methods

    
        function self = sysmat_assembler_base (Parameters)
        % Constructor.
          % Parameters: none.
        %
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            
            % nothing to be done
        end
        
        function start_assembly(self, elem_mat_nrows, elem_mat_ncols, elem_mat_nmatrices, ndofs_row, ndofs_col)
        % Start the global matrix assembly. 
        % It must be called before the first call to the method assemble.
        % elem_mat_nrows= number of rows in typical element matrix, 
        % elem_mat_ncols= number of columns in a typical element matrix, 
        % elem_mat_nmatrices= number of element matrices, 
        % ndofs_row= Total number of equations in the row direction, 
        % ndofs_col= Total number of equations in the column direction.
        %
        % No implementation in this class.
            error ('Not implemented in the abstract class');
        end

        function assemble(self, mat, dofnums_row, dofnums_col)
        % 'Assemble' a rectangular matrix.
        % The method assembles a rectangular matrix using the two vectors of 
        % degree-of-freedom numbers for the rows and columns.
        %
        % No implementation in this class.
            error ('Not implemented in the abstract class');
        end
        
        function assemble_symmetric(self, mat,  nums)
        % 'Assemble'  a square symmetric matrix.
        % The method assembles a square symmetric matrix and the vector of 
        % degree-of-freedom numbers for the rows and columns.
        %
        % No implementation in this class.
            error ('Not implemented in the abstract class');
        end
        
        function S= make_matrix (self)
        % Make a sparse matrix.
        % The method makes a sparse matrix from the assembly buffers.
        %
        % No implementation in this class.
            error ('Not implemented in the abstract class');
        end
        
    end
end

