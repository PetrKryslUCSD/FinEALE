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
        %             emask_row=dofnums_row>self.inv_dofum;
        %             emask_col=dofnums_col>self.inv_dofum;
        %             conns_row=dofnums_row(emask_row);
        %             conns_col=dofnums_col(emask_col);
        %             rmat=mat(emask_row,emask_col);
            nrows=length(dofnums_row); ncolumns=length(dofnums_col);
            ntotal =nrows*ncolumns;
            buffer_range=self.buffer_pointer:self.buffer_pointer+ntotal-1;
            self.matbuffer(buffer_range) = mat(:); % serialized matrix
            buffer_range=self.buffer_pointer:self.buffer_pointer+nrows-1;
            for   k=1:ncolumns
                self.rowbuffer(buffer_range) = dofnums_row(:);
                buffer_range=buffer_range+nrows;
            end
            buffer_range=self.buffer_pointer:self.buffer_pointer+ncolumns-1;
            for   k=1:nrows
                self.colbuffer(buffer_range) = dofnums_col(k);
                buffer_range=buffer_range+ncolumns;
            end
            self.buffer_pointer=self.buffer_pointer+ntotal;
        end
        
        function assemble_symmetric(self, mat, dofnums)
        % Assembly of a square symmetric matrix.
        % The method assembles a square symmetric matrix using the vector of 
        % equation numbers for the rows and columns.
        nrows=length(dofnums); ncolumns=nrows;
            ntotal =nrows*ncolumns;
            buffer_range=self.buffer_pointer:self.buffer_pointer+ntotal-1;
            self.matbuffer(buffer_range) = mat; % serialized matrix
            self.rowbuffer(buffer_range) = reshape (repmat(dofnums,1,ncolumns),1,ntotal);
            self.colbuffer(buffer_range) = reshape (repmat(dofnums',nrows,1),1,ntotal);
            self.buffer_pointer=self.buffer_pointer+ntotal;
        end
        
        function S= make_matrix (self)
            % Make a sparse matrix.
            % The method makes a sparse matrix from the assembly buffers.
            self.rowbuffer(self.rowbuffer(:)<1)=self.ndofs_row+1;
            self.colbuffer(self.colbuffer(:)<1)=self.ndofs_col+1;
            S = sparse(self.rowbuffer(1:self.buffer_pointer-1), ...
                self.colbuffer(1:self.buffer_pointer-1), self.matbuffer(1:self.buffer_pointer-1), ...
                self.ndofs_row+1, self.ndofs_col+1);
            S = S(1:end-1,1:end-1);% remove the extraneous rows and columns
            self.buffer_pointer=[];
            self.matbuffer=[];
            self.rowbuffer=[];
            self.colbuffer=[];
            self.buffer_pointer=[];
            self.ndofs_row= []; self.ndofs_col=  [];
        end
        
    end
end

