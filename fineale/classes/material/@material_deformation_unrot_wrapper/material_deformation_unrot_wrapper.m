classdef material_deformation_unrot_wrapper < material_deformation_triax
    %     Class that represents deformable nonlinear material in the
    %     unrotated format.
    %
    %     The wrapped material needs to be of the small-strain,
    % small-deformation kind. 
    %
    % Reference
    % WARP3D-­Release 17.6.0
    % 3-­D Dynamic Nonlinear Fracture Analyses of Solids  Using
    % Parallel Computers, Brian Healy, Arne Gullerud,  Kyle
    % Koppenhoefer, Arun Roy, Sushovan RoyChowdhury, Jason Petti, Matt 
    % Walters, Barron Bichon,  Kristine Cochran, Adam Carlyle,  James 
    % Sobotka,  Mark Messner  and Robert Dodds University	 of Illinois
    % at Urbana -­?Champaign, July 2015
    %
    
    properties
        wrapped_material=[]; % The wrapped small-strain material.
    end
    
    methods
        
        function self = material_deformation_unrot_wrapper(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            self.wrapped_material = Parameters.wrapped_material;
            self.property = self.wrapped_material.property;
        end
        
        function val = are_tangent_moduli_constant (self)
        % Is the material stiffness matrix independent of location (constant, 
        % corresponding to a homogeneous material)?
            val =false;
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
            
            % Get the basic kinematic variables
            Fn1 = context.Fn1;
            Fn  = context.Fn;
            dt  = context.dt;
            [R,U] = polardecomp(Fn1);
            
            midstep_F=compute_midstep_def_grad(Fn1,Fn);
            [midstep_R,~] = polardecomp(midstep_F);
            rel_F=compute_rel_def_grad(Fn1,Fn);
            [midstep_D,~]=compute_midstep_rate_of_def (rel_F, dt);
            
            if nargout == 1 % Pure retrieval,  no update
                context.strain =ms.strain;
                out =state(self.wrapped_material, ms.wrapped_ms, context);
            else
                % Unrotated rate of deformation
                d = midstep_R'*midstep_D*midstep_R;
                dv=self.strain_3x3t_to_6v (d);
                newms=ms;
                newms.strain = newms.strain + dt*dv;
                context.strain =newms.strain;
                [out,newwrapped_ms] =state(self.wrapped_material, ms.wrapped_ms, context);
                newms.wrapped_ms = newwrapped_ms;
            end
            
            if strcmp(context.output,'Cauchy')
                % Rotate the output stress
                sig= R*self.stress_6v_to_3x3t(out)*R';
                out=self.stress_3x3t_to_6v(sig);
            end
                      
            return;
            
        end
        
        function D = tangent_moduli(self, context)
            % Calculate the material stiffness matrix.
            %   Call as:
            %     D = tangent_moduli(m, context)
            %   where
            %     m=material
            %    context=struct
            % with fields
            %  - mandatory
            %  - optional
            %     stiff_type=type of the stiffness,
            %        either 'Lagrangean'
            %        or 'Eulerian'
            %              requires field 'F'= current deformation gradient
            % the output arguments are
            %     D=matrix 6x6
            %
            stiff_type='Eulerian';
            if isfield(context,'stiff_type')
                stiff_type= context.stiff_type;
            end
            switch (stiff_type)
                case 'Eulerian'
                    % Get the basic kinematic variables
                    Fn1 = context.Fn1;
                    Fn  = context.Fn;
                    dt  = context.dt;
                    [R,U] = polardecomp(Fn1);
                    
                    midstep_F=compute_midstep_def_grad(Fn1,Fn);
                    [midstep_R,~] = polardecomp(midstep_F);
                    rel_F=compute_rel_def_grad(Fn1,Fn);
                    [midstep_D,~]=compute_midstep_rate_of_def (rel_F, dt);
                    
                    context.ms=context.ms.wrapped_ms;
                    Dur=  tangent_moduli(self.wrapped_material, context);
                    D = rotate_stiffness(self,Dur,R');
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
        % Create a new material state.
        %
        function ms = newmatstate (self)
            ms.wrapped_ms = self.wrapped_material.newmatstate();
            ms.strain  =zeros(6,1);
        end
        
    end
    
end
