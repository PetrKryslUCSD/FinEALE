classdef material_deformation_bb_transviso_triax < material_deformation_triax
    % Class that represents triaxial transversely isotropic hyper elastic material of Bonet and Burton.
    %
    % @article{
    % author = {Bonet, J. and Burton, A. J.},
    % title = {A simple orthotropic, transversely isotropic hyperelastic constitutive equation for large strain computations},
    % journal = {Computer Methods in Applied Mechanics and Engineering},
    % volume = {162},
    % number = {1-4},
    % pages = {151-164},
    % year = {1998}
    % }
    
    properties (Access =private)
        lambda =[];
        mu =[];
        alpha =[];
        beta =[];
        gamma =[];
    end
    
    methods
        
        function self = material_deformation_bb_transviso_triax(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            %     E1 =Young's modulus in the direction of the axis of transverse isotropy (the so-called longitudinal elastic modulus)
            %     E2 =Young's modulus in the plane of transverse isotropy (the so-called transverse elastic modulus)
            %     nu12 =Poisson ratio for loading along the axis of transverse isotropy
            %     nu23 =Poisson ratio for loading in the plane orthogonal to the axis of transverse isotropy
            %     G12=Shear modulus for shear between the axis of transverse isotropy in the plane orthogonal to it
            EL = self.property.E1;
            nuLT = self.property.nu12;
            nu = self.property.nu23;
            GLT = self.property.G12;
            ET = self.property.E2;
            G =ET/(2*(1+nu));
            n =EL/ET;% Larger than one for stiff fibers
            m =1-nu-2*nu^2/n;% Note: this parameter was defined incorrectly  in BB
            self.mu=G;
            self.alpha=self.mu-GLT;
            self.beta=G/2 + (ET*nuLT)/(4*m) - (ET*(- nuLT^2 + n))/(4*m*n*(nu + 1));
            self.gamma=(ET*(- nuLT^2 + n))/(8*m*n*(nu + 1)) - (ET*nuLT)/(4*m) - (ET*n*(nu - 1))/(8*m) - GLT/2;
            self.lambda =(ET*(- nuLT^2 + n))/(m*n*(nu + 1)) - 2*G;

            % clear all
            % syms E1 E2 E3 nu12 nu13 nu23 G12 G13 G23 EL ET nuLT nu GLT G n real
            % % n=EL/ET;
            % E1=n*ET;
            % E2 =ET; E3 =ET;
            % nu12=nuLT; nu13=nuLT;
            % G12=GLT; G13=GLT;
            % nu23=nu; G23=G;
            % compliance =[1/E1      -nu12/E1    -nu13/E1  0   0   0;...
            %                 -nu12/E1     1/E2      -nu23/E2  0   0   0;...
            %                 -nu13/E1   -nu23/E2       1/E3   0   0   0;...
            %                 0           0           0 1/G12 0   0;...
            %                 0           0           0   0 1/G13 0;...
            %                 0           0           0   0   0 1/G23]
            % stiffness =simple (compliance^-1)
            % D=subs(stiffness,'(2*nuLT^2 - n + n*nu)','-m*n')
            %
            % syms lambda mu real
            %             mI = diag([1 1 1 0.5 0.5 0.5]);
            %             m1 = [1 1 1 0 0 0]';
            %           Diso=  lambda * m1 * m1' + 2 * mu * mI
            %
            % syms alpha beta gamma real
            % Dtrn=0*D; Dtrn(1,1) =8*beta+8*gamma-4*alpha;
            % Dtrn(1,2) =4*beta; Dtrn(2,1) =4*beta;
            % Dtrn(1,3) =4*beta; Dtrn(3,1) =4*beta;
            % Dtrn(4,4) =-alpha; Dtrn(5,5) =-alpha;
            %
            %           Dsst=Diso+Dtrn
            %
            % self = solve(Dsst(1,1)-D(1,1),Dsst(1,2)-D(1,2),Dsst(2,2)-D(2,2),Dsst(4,4)-D(4,4),Dsst(6,6)-D(6,6),alpha, beta, gamma,lambda, mu)
            %
            %
            % self.alpha
            %       self.beta
            %      self.gamma
            %     self.lambda
            %         self.mu
            %
            % compliance =
            %
            % [     1/(ET*n), -nuLT/(ET*n), -nuLT/(ET*n),     0,     0,   0]
            % [ -nuLT/(ET*n),         1/ET,       -nu/ET,     0,     0,   0]
            % [ -nuLT/(ET*n),       -nu/ET,         1/ET,     0,     0,   0]
            % [            0,            0,            0, 1/GLT,     0,   0]
            % [            0,            0,            0,     0, 1/GLT,   0]
            % [            0,            0,            0,     0,     0, 1/G]
            %
            % Warning: simple will be removed in a future release. Use simplify instead.
            % > In sym.simple at 41
            %
            % stiffness =
            %
            % [ (ET*n^2*(nu - 1))/(2*nuLT^2 - n + n*nu),                     -(ET*n*nuLT)/(2*nuLT^2 - n + n*nu),                     -(ET*n*nuLT)/(2*nuLT^2 - n + n*nu),   0,   0, 0]
            % [      -(ET*n*nuLT)/(2*nuLT^2 - n + n*nu),  -(ET*(- nuLT^2 + n))/((nu + 1)*(2*nuLT^2 - n + n*nu)), -(ET*(nuLT^2 + n*nu))/((nu + 1)*(2*nuLT^2 - n + n*nu)),   0,   0, 0]
            % [      -(ET*n*nuLT)/(2*nuLT^2 - n + n*nu), -(ET*(nuLT^2 + n*nu))/((nu + 1)*(2*nuLT^2 - n + n*nu)),  -(ET*(- nuLT^2 + n))/((nu + 1)*(2*nuLT^2 - n + n*nu)),   0,   0, 0]
            % [                                       0,                                                      0,                                                      0, GLT,   0, 0]
            % [                                       0,                                                      0,                                                      0,   0, GLT, 0]
            % [                                       0,                                                      0,                                                      0,   0,   0, G]
            %
            %
            % D =
            %
            % [ -(ET*n*(nu - 1))/m,                         (ET*nuLT)/m,                         (ET*nuLT)/m,   0,   0, 0]
            % [        (ET*nuLT)/m,  (ET*(- nuLT^2 + n))/(m*n*(nu + 1)), (ET*(nuLT^2 + n*nu))/(m*n*(nu + 1)),   0,   0, 0]
            % [        (ET*nuLT)/m, (ET*(nuLT^2 + n*nu))/(m*n*(nu + 1)),  (ET*(- nuLT^2 + n))/(m*n*(nu + 1)),   0,   0, 0]
            % [                  0,                                   0,                                   0, GLT,   0, 0]
            % [                  0,                                   0,                                   0,   0, GLT, 0]
            % [                  0,                                   0,                                   0,   0,   0, G]
            %
            %
            % Diso =
            %
            % [ lambda + 2*mu,        lambda,        lambda,  0,  0,  0]
            % [        lambda, lambda + 2*mu,        lambda,  0,  0,  0]
            % [        lambda,        lambda, lambda + 2*mu,  0,  0,  0]
            % [             0,             0,             0, mu,  0,  0]
            % [             0,             0,             0,  0, mu,  0]
            % [             0,             0,             0,  0,  0, mu]
            %
            %
            % Dsst =
            %
            % [ 8*beta - 4*alpha + 8*gamma + lambda + 2*mu, 4*beta + lambda, 4*beta + lambda,          0,          0,  0]
            % [                            4*beta + lambda,   lambda + 2*mu,          lambda,          0,          0,  0]
            % [                            4*beta + lambda,          lambda,   lambda + 2*mu,          0,          0,  0]
            % [                                          0,               0,               0, mu - alpha,          0,  0]
            % [                                          0,               0,               0,          0, mu - alpha,  0]
            % [                                          0,               0,               0,          0,          0, mu]
            %
            %
            % self =
            %
            %      alpha: [1x1 sym]
            %       beta: [1x1 sym]
            %      gamma: [1x1 sym]
            %     lambda: [1x1 sym]
            %         mu: [1x1 sym]
            %
            %
            % ans =
            %
            % G - GLT
            %
            %
            % ans =
            %
            % G/2 + (ET*nuLT)/(4*m) - (ET*(- nuLT^2 + n))/(4*m*n*(nu + 1))
            %
            %
            % ans =
            %
            % (ET*(- nuLT^2 + n))/(8*m*n*(nu + 1)) - (ET*nuLT)/(4*m) - (ET*n*(nu - 1))/(8*m) - GLT/2
            %
            %
            % ans =
            %
            % (ET*(- nuLT^2 + n))/(m*n*(nu + 1)) - 2*G
            %
            %
            % ans =
            %
            % G
            %
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
            
            Id=eye(3,3);
            % Finger deformation tensor
            b=Fn1*Fn1';
            J=det(Fn1);
            kirch_iso = (self.mu) * (b - eye(3,3)) + (self.lambda *log(J)) * Id;
            % The transversely  isotropic model has the first basis vector as the axis of isotropy
            a=Fn1(:,1); % this is the vector along the fiber
            I4 =a'*a; % pseudo-invariant
            kirch_trn = 2*self.beta*(I4-1)*Id ...
                + 2*(self.alpha+2*self.beta*log(J)+2*self.gamma*(I4-1))*a*a' ...
                - self.alpha*(b*a*a'+a*(b*a)');
            
            Cauchy=(kirch_iso+kirch_trn)/J;
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                switch context.output
                    case 'Cauchy'
                        out = Cauchy;
                    case 'pressure'
                        out = -(sum(diag(Cauchy))/3);
                    case 'princCauchy'
                        [V,D]=eig(Cauchy);
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        stress = self.stress_3x3t_to_6v(Cauchy);
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'equiv_pl_def'
                        out =0.0;
                    case'strain_energy'
                        out = 0.0;
                    otherwise
                        out = [];
                end
            end
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
            Fn1 = context.Fn1;
            switch (stiff_type)
                case 'Eulerian'
                    ms=context.ms;
                    D=zeros(6,6);
                    dep=sqrt(eps);
                    ix=[1,4,5,;0,2,6;0,0,3];
                    ctx.Fn1 = context.Fn1; ctx.Fn = context.Fn; 
                    ctx.dt = context.dt; ctx.dT = []; ctx.output='Cauchy';;
                    [sigmav0, ~] = state (self, ms, ctx);
                    for i=1:3
                        for j=i:3
                            Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                            ctx.Fn1 = Fn1 + Finc*Fn1;
                            [sigmav, ~] = state (self, ms, ctx);
                            D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
                        end
                    end
                    D= (D+D')/2;
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
    end
    
end
