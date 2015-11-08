classdef femm_deformation_nonlinear_h8msgso < femm_deformation_nonlinear
    % Class for large-strain nonlinear deformation based on the mean-strain 
    % hexahedra and stabilization by Gaussian quadrature.
    %
    % Reference:
    % Krysl, P.,  Mean-strain eight-node hexahedron with stabilization by
    % energy sampling, INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN
    % ENGINEERING Article first published online : 24 JUN 2014, DOI:
    % 10.1002/nme.4721
    %
    % Krysl, P.,  Optimized Energy-Sampling Stabilization of the
    % Mean-strain 8-node Hexahedron, submitted to IJNME 23 August 2014.
    %
    % Krysl, P.,  Mean-strain 8-node Hexahedron with Optimized
    % Energy-Sampling Stabilization for Large-strain Deformation, submitted
    % to IJNME 18 September 2014.
    
    properties %(Hidden, SetAccess = private)
        % Stabilization material.       
        stabilization_material  = [];
        nu=[];
        phis = [];%  stabilization parameters, computed in the update()  method
    end
    
    
    methods % constructor
        
        function self = femm_deformation_nonlinear_h8msgso (Parameters)
        % Parameters
        % stabilization_material = optional, stabilization material; if not
        %           supplied, the stabilization material is a neo-Hookean
        %           material  with material properties adjusted to avoid
        %           volumetric locking
        % nconstrained  = used by the stabilization FEMM, number of constrained modes in the material matrix
            if nargin <1
                Parameters = struct('stabilization_material',[]);
            end
            self = self@femm_deformation_nonlinear(Parameters);
            if (isfield(Parameters,'stabilization_material'))
                self.stabilization_material= Parameters.stabilization_material;
            else
                %     We will have to make our own stabilization material.
                try
                    E = self.material.property.E;% Isotropic  material?
                    if (self.material.property.nu<0.3)
                        nu=self.material.property.nu;
                    else
                        nu=0.3+ (self.material.property.nu-0.3)/2;
                    end
                catch
                    try % Orthotropic material?
                        E = self.material.property.E1;
                        try
                            E = [E,self.material.property.E2];
                        catch
                        end
                        try
                            E = [E,self.material.property.E3];
                        catch
                        end
                        E=min(E);;
                        nu=min([self.material.property.nu12,self.material.property.nu13,self.material.property.nu23]);
                    catch
                        error('Do not know how to create stabilization material');
                    end
                end
                %                 nu=0.; % Experiment in
                if (isfield( Parameters, 'match_stabilization'))
                    if (strcmp(class(self.material ),'material_deformation_stvk_triax'))
                        prop = property_deformation_linear_iso(struct('E',E,'nu',nu));
                        self.stabilization_material = material_deformation_stvk_triax(struct('property',prop));
                    else
                        prop = property_deformation_neohookean (struct('E',E,'nu',nu));
                        self.stabilization_material = material_deformation_neohookean_triax(struct('property',prop));
                    end
                else
                    prop = property_deformation_neohookean (struct('E',E,'nu',nu));
                    self.stabilization_material = material_deformation_neohookean_triax(struct('property',prop));
                end
                %
                %                 prop = property_deformation_neohookean (struct('E',E,'nu',nu));
                %                 self.stabilization_material = material_deformation_neohookean_triax(struct('property',prop));
            end
            % Now try to figure out which Poisson ratio to use in the
            % optimal scaling factor (to account for geometry)
            try % Isotropic material
                self.nu = self.material.property.nu;
            catch
                try %  Orthotropic material
                    self.nu = self.material.property.nu12;
                    try
                        self.nu = [self.nu,self.material.property.nu13];
                    catch
                    end
                    try
                        self.nu = [self.nu,self.material.property.nu23];
                    catch
                    end
                    self.nu=max(self.nu);;
                catch
                    error('Do not know how to pick Poisson ratio');
                end
            end
            %             self.nu=self.stabilization_material.property.nu; % Experiment
        end
        
    end
    
    methods (Hidden)
        
        
        function phi=stab_fraction(self,J);
            %             Compute  the stabilization parameter
            h2=diag(J'*J);
            % Jt=J';
            %             h2=[Jt(1,:)*J(:,1),Jt(2,:)*J(:,2),Jt(3,:)*J(:,3)];
            %             h2=sort(h2,'descend');
            %             h2=sort(h2,'descend');%Experiment
            Phistress= ( 2*(1+self.nu)*(min(h2)/max(h2)) );% Plane stress
            Phistrain= ( 2*(1+self.nu)*(1-self.nu)/(1+self.nu)/(1-2*self.nu)*(min(h2)/max(h2)) );%  Plane-strain
            Phi = (Phistress+Phistrain)/2;
            %             Phi = (Phistrain);
            %             Phi = (Phistress);
            phi = Phi/(1+Phi);
        end 
        
    end
    
end
