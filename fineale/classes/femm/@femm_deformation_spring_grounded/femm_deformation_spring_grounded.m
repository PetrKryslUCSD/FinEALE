classdef femm_deformation_spring_grounded < femm_base
    % Class to represent grounded spring elements.
    %
    properties
        translation_stiffness_matrix=  {};% cell array of matrices, one per finite element
    end
    
    
    methods % constructor
        
        function self = femm_deformation_spring_grounded (Parameters)
            % Constructor.
            % Parameters:
            %  translation_stiffness_matrix=  % cell array of matrices, one per finite element
            %
            % Example:
            % e1,e2,e3= ortho normal basis vectors
            % Stiffness matrix for constraints tying a node to a curve  defined by
            % e1:
            %         Gstiffness{end+1}=K*(e2*e2'+e3*e3');
%
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct('control_points',[]);
            end
            self = self@femm_base(Parameters);
            if (isfield(Parameters,'translation_stiffness_matrix'))
                self.translation_stiffness_matrix= Parameters.translation_stiffness_matrix;
            end
            if (strcmp(class(self.translation_stiffness_matrix), 'double'))
                k=self.translation_stiffness_matrix;
                self.translation_stiffness_matrix= cell(count( self.fes),1);
                for    COUNTER =1:length(self.translation_stiffness_matrix)
                    self.translation_stiffness_matrix{COUNTER}=k;
                end; clear COUNTER
            end

        end
        
    end
    
    methods
        
        
        function K = stiffness (self, assembler, geom0, u)
            % Compute the stiffness matrix for the grounded spring elements.
            %
            % function K = stiffness (self, assembler, geom0, u1, dchi)
            %
            %     assembler = descendent of the sysmat_assembler class
            %     geom0=initial geometry field; locations of nodes
            %     u1=current geometry increment field; displacements of nodes
            %     Rfield1=field of nodal frames, current deformation
            %     dchi=field of translation and rotation increments
            %             Equation numbers are taken from the field dchi.
            %
            
            conns = self.fes.conn; % connectivity
            Xs =geom0.values;
            U1s =u.values;
            X1s=Xs+U1s;
            Kedim =3*self.fes.nfens;
            start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
            for i=1:count(self.fes)                % Now loop over all fes in the block
                conn = conns(i,:); % connectivity
                dofnums =reshape(u,gather_dofnums(u,conn));
                Ke=zeros(3);
                if (~isempty(self.translation_stiffness_matrix))
                    Ke(1:3,1:3) = self.translation_stiffness_matrix{i};
                    assemble_symmetric(assembler, Ke, dofnums);
                end
            end
            K = make_matrix (assembler);
        end
        
        function [self,F] = restoring_force(self, assembler, geom0, u1, Rfield1, dchi)
            % Compute the grounded spring restoring force vector.
            %
            % function [self,F] = restoring_force(self, assembler, geom0, u1, dchi)
            %
            %     assembler = descendent of the sysvec_assembler class
            %     geom=initial geometry field; locations of nodes
            %     Rfield1=field of nodal frames, current deformation
            %     dchi=field of translation and rotation increments
            %          Equation numbers are taken from the field dchi.
            %
            conns = self.fes.conn; % connectivity
            Xs =geom0.values;
            U1s =u1.values;
            Rs = Rfield1.values;
            X1s=Xs+U1s;
            % Prepare assembly
            start_assembly(assembler, dchi.nfreedofs);
            for i=1:count(self.fes)
                conn = conns(i,:); % connectivity
                dofnums =reshape(dchi,gather_dofnums(dchi,conn));
                x0 = Xs(conn,:); % initial coordinates of nodes
                Fe=zeros(3,1);
                if (~isempty(self.translation_stiffness_matrix))
                    Ke = self.translation_stiffness_matrix{i};
                    us=U1s(conn,:)';
                    Fe(1:3)= -Ke*us;
                end
                assemble(assembler, Fe, dofnums(1:6));
            end
            F = make_vector (assembler);
        end
        
        
    end
    
end
