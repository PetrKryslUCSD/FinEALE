classdef sysvec_assembler < handle
% The class sysvec_assembler is for assembling of a system
% column vector from elementwise vectors.
%
% Note that this class is derived from the handle class.   The properties of the 
% class may be changed from the outside –- the variables of this class are passed 
% by reference not by value.
%

	properties (Constant, GetAccess = public)
	% Fixed degrees of freedom numbers are given this value:
     % it indicates that this is not a valid free  degree of freedom number.
		invalid_dofnum=sysmat_assembler_base.fixed_dofnum;
	end
    
    properties  (Hidden, GetAccess= protected, SetAccess= protected)
        F_buffer= [];
    end
    
	methods

        function self = sysvec_assembler (Parameters)
        % Constructor.
          % Parameters: none.
        %
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            
            % nothing to be done
        end
        
        function start_assembly(self, ndofs_row)
        % Start assembly.
        %
        % function start_assembly(self, ndofs_row)
        %
        % The method makes the buffer for the vector assembly. It must be called before 
        % the first call to the method assemble.
        % ndofs_row= Total number of degrees of freedom.
            self.F_buffer= zeros(ndofs_row,1)+self.invalid_dofnum;
        end

        function assemble(self, vec, dofnums)
        % Assembly of elementwise vector.
        %
        % function assemble(self, vec, dofnums)
        %
        % The method assembles a column element vector using the vector of 
        % equation numbers for the rows.
            for i = 1:length(dofnums)
                    gi = dofnums(i);
                    if (gi ~= self.invalid_dofnum)
                        self.F_buffer(gi) = self.F_buffer(gi) + vec(i);
                    end
            end
        end
        
        function F= make_vector (self)
        % Make the global vector.
        %
        % function F= make_vector (self)
        %
        % The method makes a vector from the assembly buffers.
            F=self.F_buffer;
        end
        
    end
end

