classdef sysmat_assembler_sparse_fixed_symm < sysmat_assembler_base
    % The class is for assembling of a sparse global matrix from symmetric fixed-size elementwise matrices.
    %
    
    properties (Hidden, GetAccess= protected, SetAccess= protected)
        buffer_length= [];
        elem_mat_nrows=[]; elem_mat_ncols=[];
        matbuffer=[];
        dofnumsbuffer=[];
        buffer_pointer=[];
        ndofs_row= []; ndofs_col=  [];
        inv_dofnum=sysmat_assembler_base.fixed_dofnum;
    end
    
    methods
        
        function self = sysmat_assembler_sparse_fixed_symm (Parameters)
            % Constructor.
            % Parameters: none.
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin > 0
            end
            self.buffer_length= [];
            self.elem_mat_nrows=[]; self.elem_mat_ncols=[];
            self.matbuffer=[];
            self.dofnumsbuffer=[];
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
            if ( elem_mat_nrows~=elem_mat_ncols)
                error('Square elementwise matrices assumed here')
            end
            self.buffer_length =elem_mat_nmatrices;
            self.elem_mat_nrows=elem_mat_nrows; self.elem_mat_ncols=elem_mat_ncols;
            self.dofnumsbuffer = zeros (elem_mat_ncols,self.buffer_length);
            self.matbuffer = zeros (elem_mat_nrows*elem_mat_ncols,self.buffer_length);
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
            self.matbuffer(:,self.buffer_pointer) = mat(:); % serialized matrix
            self.dofnumsbuffer(:,self.buffer_pointer) = dofnums_col;
            self.buffer_pointer=self.buffer_pointer+1;
        end
        
        function assemble_symmetric(self, mat, dofnums)
            % Assembly of a square symmetric matrix.
            % The method assembles a square symmetric matrix using the vector of
            % equation numbers for the rows and columns.
            self.matbuffer(:,self.buffer_pointer) = mat(:); % serialized matrix
            self.dofnumsbuffer(:,self.buffer_pointer) = dofnums;
            self.buffer_pointer=self.buffer_pointer+1;
        end
        
        function S= make_matrix (self)
            % Make a sparse matrix.
            % The method makes a sparse matrix from the assembly buffers.
            self.buffer_pointer=self.buffer_pointer-1;
            self.dofnumsbuffer=self.dofnumsbuffer(:,1:self.buffer_pointer);
            self.matbuffer=self.matbuffer(:,1:self.buffer_pointer);
            self.dofnumsbuffer(self.dofnumsbuffer(:)<1)=self.ndofs_col+1;
            B = repmat(self.dofnumsbuffer',1,1,size(self.dofnumsbuffer,1));
            rb = permute(B,[2,3,1]);
            cb = permute(B,[3,2,1]);
            clear B
            S = sparse(rb(:), cb(:), self.matbuffer(:), self.ndofs_row+1, self.ndofs_col+1);
            S = S(1:end-1,1:end-1);% remove the extraneous rows and columns
            self.buffer_pointer=[];
            self.elem_mat_nrows=[]; self.elem_mat_ncols=[];
            self.matbuffer=[];
            self.dofnumsbuffer=[];
            self.buffer_pointer=[];
            self.ndofs_row= []; self.ndofs_col=  [];
        end
        
    end
end

