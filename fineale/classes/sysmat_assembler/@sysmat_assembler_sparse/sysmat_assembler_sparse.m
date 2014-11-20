classdef sysmat_assembler_sparse < sysmat_assembler_base
% The class is for assembling of a sparse global matrix from elementwise matrices.
% 

    properties (Hidden, GetAccess= protected, SetAccess= protected)
        buffer_length= [];
        matbuffer=[];
        rowbuffer=[];
        colbuffer=[];
        buffer_pointer=[];
        ndofs_row= []; ndofs_col=  [];
        inv_dofnum=sysmat_assembler_base.fixed_dofnum;
    end
    
	methods
        
        function self = sysmat_assembler_sparse (Parameters)
        % Constructor.
          % Parameters: none.
        %
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin > 0
            end
            self.buffer_length= [];
            self.matbuffer=[];
            self.rowbuffer=[];
            self.colbuffer=[];
            self.buffer_pointer=[];
            self.ndofs_row= []; 
            self.ndofs_col=  [];
        end
        
        function start_assembly(self, elem_mat_nrows, elem_mat_ncols, elem_mat_nmatrices, ndofs_row, ndofs_col)
        % Start the assembly of a global matrix.
        % The method makes buffers for matrix assembly. It must be called before 
        % the first call to the method assemble.
        % elem_mat_nrows= number of rows in typical element matrix, 
        % elem_mat_ncols= number of columns in a typical element matrix, 
        % elem_mat_nmatrices= number of element matrices, 
        % ndofs_row= Total number of equations in the row direction, 
        % ndofs_col= Total number of equations in the column direction.
            self.buffer_length =elem_mat_nmatrices*elem_mat_nrows*elem_mat_ncols; 
            self.rowbuffer = zeros (1,self.buffer_length);
            self.colbuffer = zeros (1,self.buffer_length);
            self.matbuffer = zeros (1,self.buffer_length);
            self.buffer_pointer = 1;
            self.ndofs_row =ndofs_row;
            self.ndofs_col =ndofs_col;
        end

        function assemble(self, mat, dofnums_row, dofnums_col)
        % Assembly of a rectangular matrix.
        % The method assembles a rectangular matrix using the two vectors of 
        % equation numbers for the rows and columns.
            emask_row=dofnums_row>self.inv_dofum;
            emask_col=dofnums_col>self.inv_dofum;
            conns_row=dofnums_row(emask_row);
            conns_col=dofnums_col(emask_col);
            rmat=mat(emask_row,emask_col);
            nrows=length(conns_row); ncolumns=length(conns_col);
            ntotal =nrows*ncolumns;
            buffer_range=self.buffer_pointer:self.buffer_pointer+ntotal-1;
            self.matbuffer(buffer_range) = reshape (rmat,1,ntotal); % serialized matrix
            self.rowbuffer(buffer_range) = reshape (repmat(conns_row,1,ncolumns),1,ntotal);
            self.colbuffer(buffer_range) = reshape (repmat(conns_col',nrows,1),1,ntotal);
            self.buffer_pointer=self.buffer_pointer+ntotal;
        end
        
        function assemble_symmetric(self, mat, dofnums)
        % Assembly of a square symmetric matrix.
        % The method assembles a square symmetric matrix using the vector of 
        % equation numbers for the rows and columns.
            emask_row=dofnums >self.inv_dofnum;
            conns_row=dofnums(emask_row);
            rmat=mat(emask_row,emask_row);
            nrows=length(conns_row); ncolumns=nrows;
            ntotal =nrows*nrows;
            tempinds = conns_row(:, ones(ncolumns, 1));
            buffer_range=self.buffer_pointer:self.buffer_pointer+ntotal-1;
            self.matbuffer(buffer_range) = reshape (rmat,1,ntotal);%rmat; %
            self.rowbuffer(buffer_range) = reshape(tempinds', [], 1);
            self.colbuffer(buffer_range) = reshape(tempinds, [], 1);
            self.buffer_pointer=self.buffer_pointer+ntotal;
        end
        
        function S= make_matrix (self)
        % Make a sparse matrix.
        % The method makes a sparse matrix from the assembly buffers.
            S = sparse(self.rowbuffer(1:self.buffer_pointer-1), ...
                self.colbuffer(1:self.buffer_pointer-1), self.matbuffer(1:self.buffer_pointer-1), ...
                self.ndofs_row, self.ndofs_col);
            self.buffer_pointer=[];
            self.matbuffer=[];
            self.rowbuffer=[];
            self.colbuffer=[];
            self.buffer_pointer=[];
            self.ndofs_row= []; self.ndofs_col=  [];    
        end
        
    end
end

