classdef material_deformation < material_base
    % Class that represents deformable materials.
    %
    % The purpose of the class is to provide a variety of conversion methods for
    % stress and strain to all descendents of this class.
    
    properties
        % Mechanical property object
        % Mechanical property object which must support the following methods:
        %     are_tangent_moduli_constant();
        %     tangent_moduli();
        %     thermal_expansion_coefficients();
        %     update().
        property = [];
    end
    
    methods
        
        function self = material_deformation(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_base(Parameters);
            if isfield(Parameters,'property')
                self.property = Parameters.property;
            end
        end
        
        function v = strain_2x2t_3v (self,t)
            % Convert a matrix of 2x2 strain components  into a 3-component vector.
            %
            % function v = strain_2x2t_3v (self,t)
            %
            % Convert a matrix of 2x2 strain components (tensor)
            % into a 3-component vector.
            %
            v =zeros(3,1);
            v(1,1) = t(1,1);
            v(2,1) = t(2,2);
            v(3,1) = t(1,2) + t(2,1);
        end
        
        function t = strain_3v_to_2x2t (self,v)
            % Convert a strain 3-vector to a  matrix of 2x2 strain components (tensor)
            %
            % function t = strain_3v_to_2x2t (self,v)
            %
            % Convert a strain 3-vector to a *symmetric*
            % matrix of 2x2 strain components (tensor)
            %
            t =zeros(2,2);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(1,2) = v(3)/2;
            t(2,1) = v(3)/2;
        end
        
        function v = strain_3x3t_to_6v (self,t)
            % Convert a matrix of 3x3 strain components to a 6-component vector.
            %
            % function v = strain_3x3t_to_6v (self,t)
            %
            % convert a matrix of 3x3 strain components (tensor)
            % into a 6-component vector.
            %
            v=zeros(6,1);
            v(1,1) = t(1,1);
            v(2,1) = t(2,2);
            v(3,1) = t(3,3);
            v(4,1) = t(1,2) + t(2,1);
            v(5,1) = t(1,3) + t(3,1);
            v(6,1) = t(3,2) + t(2,3);
        end
        
        function t = strain_6v_to_3x3t (self,v)
            % Convert a strain 6-vector to a  matrix of 3x3 strain components (tensor)
            %
            % function t = strain_6v_to_3x3t (self,v)
            %
            % convert a strain 6-vector to a *symmetric*
            % matrix of 3x3 strain components (tensor)
            %
            t =zeros(3,3);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(3,3) = v(3);
            t(1,2) = v(4)/2;
            t(2,1) = v(4)/2;
            t(1,3) = v(5)/2;
            t(3,1) = v(5)/2;
            t(3,2) = v(6)/2;
            t(2,3) = v(6)/2;
        end
        
        
        function v = stress_2x2_to_3v (self,t)
            % Convert a symmetric matrix of 2x2 stress components to a 3-component vector.
            %
            % function v = stress_2x2_to_3v (self,t)
            %
            % Convert a symmetric matrix of 2x2 stress components (tensor)
            % into a 3-component vector.
            %
            v =zeros(3,1);
            v(1,1) = t(1,1);
            v(2,1) = t(2,2);
            v(3,1) = 1/2*(t(1,2) + t(2,1));
        end
        
        function t = stress_3v_to_2x2t (self,v)
            % Convert a 3-vector to a  matrix of 2x2 stress components (tensor)
            %
            % function t = stress_3v_to_2x2t (self,v)
            %
            % Convert a 3-vector to a *symmetric*
            % matrix of 2x2 stress components (tensor)
            %
            t =zeros(2,2);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(1,2) = v(3);
            t(2,1) = v(3);
        end
        
        function t = stress_3v_to_3x3t (self,v)
            % Convert a 3-vector to a matrix of 3x3 stress components (tensor)
            %
            % function t = stress_3v_to_3x3t (self,v)
            %
            % Convert a 3-vector to a *symmetric*
            % matrix of 3x3 stress components (tensor)
            %
            t=zeros (3, 3);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(1,2) = v(3);
            t(2,1) = v(3);
        end
        
        function t = stress_4v_to_3x3t (self,v)
            % Convert a 4-vector to a  matrix of 3x3 stress components (tensor).
            %
            % function t = stress_4v_to_3x3t (self,v)
            %
            %
            % Convert a 4-vector to a *symmetric*
            % matrix of 3x3 stress components (tensor).  This is
            % conversion routine that would be useful for plane strain or 
            % axially symmetric conditions.
            % The stress vector components need to be ordered as:
            %     sigma_x, sigma_y, tau_xy, sigma_z,
            % which is the ordering used for the plane-strain model reduction.
            % Therefore, for axially symmetric analysis the components need to be
            % reordered, as from the constitutive equation they come out 
            % as sigma_x, sigma_y, sigma_z, tau_xy.
            %
            t=zeros (3, 3);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(1,2) = v(3);
            t(2,1) = v(3);
            t(3,3) = v(4);
        end
        
        function t = stress_6v_to_3x3t (self,v)
            % Convert a 6-vector to a  matrix of 3x3 stress components (tensor)
            %
            % function t = stress_6v_to_3x3t (self,v)
            %
            % convert a 6-vector to a *symmetric*
            % matrix of 3x3 stress components (tensor)
            %
            t=zeros(3,3);
            t(1,1) = v(1);
            t(2,2) = v(2);
            t(3,3) = v(3);
            t(1,2) = v(4);
            t(2,1) = v(4);
            t(1,3) = v(5);
            t(3,1) = v(5);
            t(3,2) = v(6);
            t(2,3) = v(6);
        end
        
        
        function v = stress_3x3t_to_6v (self,t)
            % Convert a matrix of 3x3 stress components to a 6-component vector.
            %
            % function v = stress_3x3t_to_6v (self,t)
            %
            % Convert a matrix of 3x3 stress components (tensor)
            % into a 6-component vector.
            %
            %             v=zeros(6,1);
            %             v(1,1) = t(1,1);
            %             v(2,1) = t(2,2);
            %             v(3,1) = t(3,3);
            %             v(4,1) = 1/2*(t(1,2) + t(2,1));
            %             v(5,1) = 1/2*(t(1,3) + t(3,1));
            %             v(6,1) = 1/2*(t(3,2) + t(2,3));
            v = t([1,5,9,2,3,6])';
        end
        
        function t = strain_9v_to_6v (self,v)
            % Convert a strain 9-vector to a  strain 6-vector components (tensor)
            %
            % function t = strain_9v_to_6v (self,v)
            %
            %
            t=zeros(6,1);
            t(1,1) = v(1);
            t(2,1) = v(2);
            t(3,1) = v(3);
            t(4,1) = v(4)+v(5);
            t(5,1) = v(8)+v(9);
            t(6,1) = v(6)+v(7);
        end
        
        function t = strain_6v_to_9v (self,v)
            % Convert a strain 6-vector to a  strain 9-vector components (tensor)
            %
            % function t = strain_6v_to_9v (self,v)
            %
            %
            t=zeros(9,1);
            t(1) = v(1);
            t(2) = v(2);
            t(3) = v(3);
            t(4) = v(4)/2;
            t(5) = v(4)/2;
            t(6) = v(6)/2;
            t(7) = v(6)/2;
            t(8) = v(5)/2;
            t(9) = v(5)/2;
        end
        
        function t = stress_9v_to_6v (self,v)
            % Convert a strain 9-vector to a  strain 6-vector components (tensor)
            %
            % function t = strain_9v_to_6v (self,v)
            %
            %
            t=zeros(6,1);
            t(1,1) = v(1);
            t(2,1) = v(2);
            t(3,1) = v(3);
            t(4,1) = v(4);
            t(5,1) = v(8);
            t(6,1) = v(6);
        end
        
        function v = tensor_3x3t_double_contraction(A,B)
            % Compute a tensor double contraction (scaler).
            %
            % function v = tensor_3x3t_double_contraction(A,B)

            v = sum(sum(A.*B));
        end
         
        function cout=Lagrangean_to_Eulerian(self, C, F)
            % Convert a Lagrangean constitutive matrix to an Eulerian one.
            %
            % function cout=Lagrangean_to_Eulerian(self, C, F)
            %
            % Convert a Lagrangean constitutive matrix to an Eulerian one. NOTE:
            % the Lagrangean matrix is presumed symmetric.
            % self = object
            % C    = Lagrangean constitutive matrix, 6x6, symmetric
            % F    = current deformation gradient, F_ij = \partial x_i / \partial X_j
            error('Base class does not implement this behaviour!');
        end

        function val = are_tangent_moduli_constant (self)
            % Is the material stiffness matrix independent of location?
            % Is the material stiffness matrix independent of location (constant,
            % corresponding to a homogeneous material)?
            % No implementation in the base class.
            error('Base class does not implement this behaviour!');
        end
        
        function val = tangent_moduli(self, context)
            % Calculate the material stiffness matrix.
            %
            % function val = tangent_moduli(self, context)
            %
            %   Call as:
            %     D = tangent_moduli(m, context)
            %   where
            %     m=material
            %    context=structure
            %  with mandatory and optional fields that are required by the specific material implementation.
            %
            % the output arguments are
            %     D=matrix 6x6 in the material coordinates system
            %
            %   This is the base class: descendant needs to override this method.
            %
            % No implementation in the base class.
            error('Base class does not implement this behaviour!');
        end
        
        function [out, newms] = state(self, ms, context)
            % Retrieve material state. 
            %
            %   function [out, newms] = state(self, ms, context)
            %
            % The method is used with either one output argument or with
            % two output arguments.
            % 1.  With only "out" as output argument:  Retrieve the
            %     the requested variable from the material state. 
            % 2.  With both output arguments: Update the material state,
            %     and return the requested variable from the material state 
            %     and the updated material state.
            %
            %     The requested quantity may or may not be supported by
            %     the particular material model.(default is the stress). 
            %
            %   Input arguments:
            % self=material
            % ms = material state
            % context=structure
            %    with mandatory fields
            %       Fn1= current deformation gradient (at time t_n+1)
            %       Fn=  previous converged deformation gradient (at time t_n)
            %       dt= time step, t_n+1-t_n
            %    and optional fields
            %       output=type of quantity to output, and interpreted by the
            %           particular material; [] is returned when the material 
            %           does not recognize the requested quantity to indicate 
            %           uninitialized value.  It can be tested with isempty().
            %           output ='Cauchy' - Cauchy stress; this is the default
            %              when output type is not specified.
            %           output ='2ndPK' - 2nd Piola-Kirchhoff stress;
            %                  
            %              output ='strain_energy'
            %    It is assumed that stress is output in 6-component vector
            %    form. 
            %
            %   Output arguments:
            % out=requested quantity
            % newms=new material state; don't forget that if the update is
            %       final the material state newms must be assigned and
            %       stored.  Otherwise the material update is lost!
            %
            
            % This is the case were the material state needs to be updated
            % and then the variable requested needs to be returned.
        error('Base class does not implement this behaviour!');
        end
        
        function v = newmatstate (self)
        % Create a new material state.
        %
        % function v = newmatstate (self)
        %
        %   Call as:
        %     v = newmatstate(m)
        %   where
        %     m=an instance of a descendant of the class material
        %
        % The default is to assume a memoryless material– there is no material state.
        
            v= [];%  no data stored for memoryless material
        end
        
    end
    
end
