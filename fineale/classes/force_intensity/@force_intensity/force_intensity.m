classdef force_intensity
    % Distributed force (force intensity) class.
    %
    % The force intensity class. The physical units are
    % force per unit volume, where volume depends on to which manifold
    % the force is applied:
    % force/length^3 (when applied to a 3-D solid),
    % force/length^2 (when applied to a surface),
    % force/length^1 (when applied along a curve),
    % or force/length^0 (when applied at a point).
    %
    %
    
    properties
        % Magnitude of the distributed force.
        % Either a constant vector or a function handle.
        magn = [];
    end
    
    properties (Hidden,  GetAccess = protected, SetAccess = private)
        is_uniform  = true;
        expect_J= true;
    end
    
    methods
        
        function self = force_intensity (Parameters)
            % Constructor.
            % Parameters:
            %      magn=vector of magnitudes, dimension must correspond to the
            %            dimension of the displacement/geometry fields; magn may be an array
            %            of doubles (representing a constant, or uniform, force intensity)
            % or a function handle/in-line function with signature
            %                     function val=f(xyz)
            %            or
            %                     function val=f(xyz,J)
            %            where
            %                     xyz = location, at which the force intensity
            %                           is to be evaluated.
            %                     J = Jacobian matrix at the location above.
            %                        The Jacobian matrix could be useful for instance
            %                         for the calculation of the normal to the surface.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            self.magn = Parameters.magn;
            self.expect_J= false;
            if strcmp(class (self.magn),'inline') || strcmp(class (self.magn),'function_handle')
                for j=1:3
                    try
                        val=feval(self.magn,zeros(1,j),eye(j));
                        self.expect_J= true;
                        break;
                    catch
                        self.expect_J= ~true;
                    end
                end
            end
            if strcmp(class (self.magn),'inline') || strcmp(class (self.magn),'function_handle')
                self.is_uniform = false;
            end
        end
        
        function val = get_magn(self,xyz,J)
            % Get force intensity magnitude at the specified point.
            %
            % function val = get_magn(self,xyz,J)
            %
            %  xyz=the location at which the force intensity needs to be evaluated
            %  J= optional, Jacobian matrix of the geometric cell at which the force
            %      intensity is evaluated. This could be useful for instance
            %      for the calculation of the normal to the surface.
            if self.is_uniform
                val = self.magn;%  position- and orientation independent
            else
                if self.expect_J
                    val=feval(self.magn,xyz,J);
                else
                    val=feval(self.magn,xyz);
                end
                if length(val)==3
                    val=reshape(val, 3, 1);
                end
            end
        end
        
    end
    
end

