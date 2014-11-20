classdef femm_deformation < femm_base
    % Class for the deformation model.
    %
    
    properties
        N= [];% Normal to the FE set (meaningful only for manifold dimension <3)
        % Coefficient for surface normal spring boundary condition, in units of
         % force per unit area per unit length
        surface_normal_spring_coefficient= []; 
    end
    
    properties (Hidden, SetAccess = protected)
        hBlmat= [];% handle to the strain-displacement matrix function
    end
    
    
    methods % constructor
        
        function self = femm_deformation (Parameters)
            % Constructor.
            % Parameters: those recognized by model_base plus
            %    N = Outer normal, either as an array of numbers, or as a function handle.
            %     Signature as in this example (tangent is a matrix with tangent
            %     vectors as columns).
            %     function n= ABC_normal(x,tangent)
            %         n= cross(tangent(:,1),tangent(:,2));
            %         n=n/norm(n);
            %     end
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_base(Parameters);
            self.hBlmat= [];
            self.N= [];%
            if (nargin <1)  || isempty(Parameters)
                return
            end
            if isfield(Parameters,'N')
                self.N = Parameters.N;
            end
            if isfield(Parameters,'surface_normal_spring_coefficient')
                self.surface_normal_spring_coefficient = Parameters.surface_normal_spring_coefficient;
            end
            % Decide on which strain-displacement matrix is applicable.
            % The manifold dimension and possibly axial symmetry need to be considered.
            if (isempty(self.fes))
                error('Need at least one finite element!');
            end
            switch self.fes.dim
                case 0
                    self.hBlmat=[];
                case 1
                    self.hBlmat=@Blmat1;
                    if (isempty(self.Rm)),self.Rm=eye(1);end % for identity transformation
                case 2
                    if self.fes.axisymm
                        self.hBlmat=@Blmat2axisymm;
                    else
                        self.hBlmat=@Blmat2;
                    end
                    if (isempty(self.Rm)),self.Rm=eye(2);end % for identity transformation
                case 3
                    self.hBlmat=@Blmat3;
                    if (isempty(self.Rm)),self.Rm=eye(3);end % for identity transformation
                otherwise
                    error (' Not implemented');
            end
        end
        
    end
    
    methods (Hidden, Access = protected)
        
        function B = Blmat1(self,N,Ndersp,c,Rm)
            % Compute the strain-displacement matrix for a one-manifold element.
            %
            % B = Blmat1(self,N,Ndersp,c,Rm)
            %
            % Compute the linear, displacement independent, strain-displacement matrix
            % for a one-manifold element.
            %
            % Arguments
            %      self=finite element model
            %      N= matrix of basis function values
            %      Ndersp=matrix of basis function gradients with respect to the
            %         Cartesian coordinates in the directions of the material orientation
            %      c=array of spatial coordinates of the evaluation point
            %         in the global Cartesian coordinates
            %      Rm=orthogonal matrix with the unit basis vectors of the local
            %         material orientation coordinate system as columns;
            %         supply Rm as empty when there is no need for global-to-local
            %         transformation (for instance for isotropic materials).
            %         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
            %         of the embedding space, greater than or equal to  2.
            % Output:
            %      B = strain-displacement matrix, size (B) = [1,nfens*ndim],
            %         where  nfens = number of
            %         finite element nodes on the element.
            %
            nfn= size(Ndersp,1);
            dim =size(c,2);
            B = zeros(1,nfn*dim);
            if (isempty(Rm))% there is no global-to-local transformation
                for i= 1:nfn
                    B(:,dim*(i-1)+1:dim*i)=  Ndersp(i,1);
                end
            else% global-to-local transformation is requested
                for i= 1:nfn
                    B(:,dim*(i-1)+1:dim*i)=  Ndersp(i,1) *Rm(:,1)' ;
                end
            end
        end
        
        function B = Blmat2(self, N,Ndersp,c,Rm)
            % Compute the strain-displacement matrix for a two-manifold element.
            %
            % B = Blmat2(self,N,Ndersp,c,Rm)
            %
            % Compute the linear, displacement independent, strain-displacement matrix
            % for a two-manifold element.
            %
            % Arguments
            %      self=finite element model
            %      N= matrix of basis function values
            %      Ndersp=matrix of basis function gradients with respect to the
            %         Cartesian coordinates in the directions of the material orientation
            %      c=array of spatial coordinates of the evaluation point
            %         in the global Cartesian coordinates
            %      Rm=orthogonal matrix with the unit basis vectors of the local
            %         material orientation coordinate system as columns;
            %         supply Rm as empty when there is no need for global-to-local
            %         transformation (for instance for isotropic materials).
            %         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
            %         of the embedding space, greater than or equal to  2.
            % Output:
            %      B = strain-displacement matrix, size (B) = [3,nfens*2],
            %         where  nfens = number of
            %         finite element nodes on the element.
            %
            nfn= size(Ndersp,1);
            dim =size(c,2);
            B = zeros(3,nfn*dim);
            if (isempty(Rm))% there is no global-to-local transformation
                % The structure of the strain-displacement matrix is as
                % shown in the  snippet of the code below
                %                 for i= 1:nfn
                %                     B(:,dim*(i-1)+1:dim*i)=...
                %                         [Ndersp(i,1) 0; ...
                %                         0           Ndersp(i,2); ...
                %                         Ndersp(i,2) Ndersp(i,1) ];
                %                 end
                % This code is twice as fast
                B(1,1:dim:end) =Ndersp(:,1)';
                B(2,2:dim:end) =Ndersp(:,2)';
                B(3,1:dim:end) =B(2,2:dim:end);
                B(3,2:dim:end) =B(1,1:dim:end);
             else% global-to-local transformation is requested
                RmT=Rm(:,1:2)';
                for i= 1:nfn
                    B(:,dim*(i-1)+1:dim*i)=...
                        [Ndersp(i,1) 0; ...
                        0           Ndersp(i,2); ...
                        Ndersp(i,2) Ndersp(i,1) ]*RmT;
                end
            end
            return;
        end
        
        function B = Blmat2axisymm(self,N,Ndersp,c,Rm)
            % Compute the strain-displacement matrix for a two-manifold element for
            % axially symmetric models.
            %
            % B = Blmat2axisymm(self,N,Ndersp,c,Rm)
            %
            % Compute the linear, displacement independent, strain-displacement matrix
            % for a two-manifold element for axially symmetric models.
            %
            % Arguments
            %      self=finite element model
            %      N= matrix of basis function values
            %      Ndersp=matrix of basis function gradients with respect to the
            %         Cartesian coordinates in the directions of the material orientation
            %      c=array of spatial coordinates of the evaluation point
            %         in the global Cartesian coordinates
            %      Rm=orthogonal matrix with the unit basis vectors of the local
            %         material orientation coordinate system as columns;
            %         supply Rm as empty when there is no need for global-to-local
            %         transformation (for instance for isotropic materials).
            %         size(Rm)= [ndim,2], where ndim = number of spatial dimensions
            %         of the embedding space, greater than or equal to  2.
            % Output:
            %      B = strain-displacement matrix, size (B) = [4,nfens*2],
            %         where  nfens = number of
            %         finite element nodes on the element.
            %
            nfn= size(Ndersp,1);
            r=c(1); if r==0,r=eps; end
            dim =size(c,2);
            B = zeros(4,nfn*dim);
            if (isempty(Rm))% there is no global-to-local transformation
                for i= 1:nfn
                    B(:,dim*(i-1)+1:dim*i)=...
                        [Ndersp(i,1) 0; ...
                        0           Ndersp(i,2); ...
                        N(i)/r 0; ...
                        Ndersp(i,2) Ndersp(i,1) ];
                end
            else% global-to-local transformation is requested
                RmT=Rm(:,1:2)';
                for i= 1:nfn
                    B(:,dim*(i-1)+1:dim*i)=...
                        [Ndersp(i,1) 0; ...
                        0           Ndersp(i,2); ...
                        N(i)/r 0; ...
                        Ndersp(i,2) Ndersp(i,1) ]*RmT;
                end
            end
        end
        
        function B = Blmat3(self,N,Ndersp,c,Rm)
            % Compute the strain-displacement matrix for a three-manifold element.
            %
            % B = Blmat3(self,N,Ndersp,c,Rm)
            %
            % Compute the linear, displacement independent, strain-displacement matrix
            % for a three-manifold element.
            %
            % Arguments
            %      self=finite element model
            %      N= matrix of basis function values
            %      Ndersp=matrix of basis function gradients with respect to the
            %         Cartesian coordinates in the directions of the material orientation
            %      c=array of spatial coordinates of the evaluation point
            %         in the global Cartesian coordinates
            %      Rm=orthogonal matrix with the unit basis vectors of the local
            %         material orientation coordinate system as columns;
            %         supply Rm as empty when there is no need for global-to-local
            %         transformation (for instance for isotropic materials).
            %         size(Rm)= [ndim,ndim], where ndim = number of spatial dimensions
            %         of the embedding space.
            % Output:
            %      B = strain-displacement matrix, size (B) = [nstrain,nfens*dim],
            %         where nstrain= number of strains, dim = Number of spatial
            %         dimensions of the embedding space, and nfens = number of
            %         finite element nodes on the element.
            %
            nfn= size(Ndersp,1);
            B = zeros(6,nfn*3); %initialize
            if (isempty(Rm)) % there is no global-to-local transformation
                for i= 1:nfn
                    k=3*(i-1);
                    B(1,k+1)= Ndersp(i,1);
                    B(2,k+2)= Ndersp(i,2);
                    B(3,k+3)= Ndersp(i,3) ;
                    B(4,k+1)= Ndersp(i,2); B(4,k+2)= Ndersp(i,1);
                    B(5,k+1)= Ndersp(i,3); B(5,k+3)= Ndersp(i,1);
                    B(6,k+2)= Ndersp(i,3); B(6,k+3)= Ndersp(i,2);
                end
            else % global-to-local transformation is requested
                RmT=Rm';
                for i= 1:nfn
                    B(:,3*(i-1)+1:3*i)...
                        = [ Ndersp(i,1) 0           0  ; ...
                        0           Ndersp(i,2) 0  ; ...
                        0           0           Ndersp(i,3) ; ...
                        Ndersp(i,2) Ndersp(i,1) 0  ; ...
                        Ndersp(i,3) 0           Ndersp(i,1) ; ...
                        0           Ndersp(i,3) Ndersp(i,2) ]*RmT;
                end
            end
            return;
        end
        
        function B = Blmat3_dilat(self,N,Ndersp,c,Rm)
            % Compute the dilatation strain-displacement matrix for a three-manifold element.
            %
            % B = Blmat3_dilat(self,N,Ndersp,c,Rm)
            %
            % Compute the linear, displacement independent, dilatation strain-displacement matrix
            % for a three-manifold element.
            %
            %   Arguments
            %      self=Model
            %      N= matrix of basis function values
            %      Ndersp=matrix of basis function gradients
            %      c=array of spatial coordinates of the evaluation point
            %      Rm=orthogonal matrix with the unit basis vectors of the local
            %         coordinate system as columns; supply Rm as empty when there
            %         is no need for global-to-local transformation (for instance for
            %         isotropic materials).
            %         size(Rm)= [ndim,ndim], where ndim = Number of spatial dimensions
            %         of the embedding space.
            % Output: B = strain-displacement matrix, size (B) = [nstrain,nfens*dim],
            %      where nstrain= number of strains, dim = Number of spatial dimensions of
            %      the embedding space, and nfens = number of finite element nodes on
            %      the element.
            %
            nfn= size(Ndersp,1);
            B = zeros(6,nfn*3); %initialize
            if (isempty(Rm)) % there is no global-to-local transformation
                for i= 1:nfn
                    k=3*(i-1);
                    B(1,k+1)= Ndersp(i,1);
                    B(2,k+2)= Ndersp(i,2);
                    B(3,k+3)= Ndersp(i,3) ;
                    B(4,k+1)= Ndersp(i,2); B(4,k+2)= Ndersp(i,1);
                    B(5,k+1)= Ndersp(i,3); B(5,k+3)= Ndersp(i,1);
                    B(6,k+2)= Ndersp(i,3); B(6,k+3)= Ndersp(i,2);
                end
            else % global-to-local transformation is requested
                RmT=Rm';
                for i= 1:nfn
                    B(:,3*(i-1)+1:3*i)...
                        = [ Ndersp(i,1) Ndersp(i,2) Ndersp(i,3); ...
                        Ndersp(i,1) Ndersp(i,2) Ndersp(i,3); ...
                        Ndersp(i,1) Ndersp(i,2) Ndersp(i,3); ...
                        0           0           0  ; ...
                        0           0           0  ; ...
                        0           0           0  ]*RmT;
                end
            end
            B =B/3;
            return;
        end
        
        
        function norml = normal(self,c,J)
            % Compute local normal. This makes sense for bounding surfaces only.
            %
            % function norml = normal(self,c,J)
            %
            %      self=finite element block
            %      c=spatial location
            %      J= Jacobian matrix
            %
            norml=self.N;
            if isempty(norml)
                % Produce a default normal
                if (size(J,1)==3)&&(size(J,2)==2)% surface in three dimensions
                    norml= skewmat(J(:,1))*J(:,2);% outer normal to the surface
                    norml=norml/norm(norml);
                elseif (size(J,1)==2)  && (size(J,2)==1)% curve in two dimensions
                    norml= [J(2,1);-J(1,1)];% outer normal to the contour
                    norml=norml/norm(norml);
                else
                    error(['No definition of normal vector']);                
                end
            else
                if  strcmp(class(norml),'function_handle')
                    norml= feval(norml, c, J);
                end
            end
        end
        
        
    end
    
end



% function B = divmat1(self,Ndersp,c)
% nfn= size(Ndersp,1);
% dim =size(Ndersp,2);
% B = zeros(1,nfn*dim);
% if (isempty(Rm))% there is no global-to-local transformation
% for i= 1:nfn
% B(:,dim*(i-1)+1:dim*i)=  Ndersp(i,1);
% end
% else% global-to-local transformation is requested
% for i= 1:nfn
% B(:,dim*(i-1)+1:dim*i)=  Ndersp(i,1) *Rm(:,1)' ;
% end
% end
% end

% function B = divmat2(self, N,Ndersp,c)
% nfn= size(Ndersp,1);
% dim =size(Ndersp,2);
% B = zeros(1,nfn*dim);
% if (isempty(Rm))% there is no global-to-local transformation
% for i= 1:nfn
% B(:,2*(i-1)+1:2*i)= [Ndersp(i,1)    Ndersp(i,2) ];
% end
% else% global-to-local transformation is requested
% RmT=Rm(:,1:2)';
% for i= 1:nfn
% B(:,2*(i-1)+1:2*i)= [Ndersp(i,1)    Ndersp(i,2) ] *RmT ;
% end
% end
% end

% function B = divmat2axisymm(self,Ndersp,c)
% nfn= size(Ndersp,1);
% r=c(1); if r==0,r=eps; end
% B = zeros(2,nfn*dim);
% for i= 1:nfn
% B(:,dim*(i-1)+1:dim*i)=...
% [Ndersp(i,1) 0; ...
% 0           Ndersp(i,2); ...
% N(i)/r 0]*Rm(:,1:2)';
% end
% return;
% end

% function B = divmat3(self,N,Ndersp,c)
% nfn= size(Ndersp,1);
% B = zeros(1,nfn*3);
% if (isempty(Rm)) % there is no global-to-local transformation
% for i= 1:nfn
% B(:,3*(i-1)+1:3*i)  = [ Ndersp(i,1)  Ndersp(i,2)  Ndersp(i,3) ];
% end
% else % global-to-local transformation is requested
% RmT=Rm';
% for i= 1:nfn
% B(:,3*(i-1)+1:3*i)  = [ Ndersp(i,1)  Ndersp(i,2)  Ndersp(i,3) ]*RmT;
% end
% end
% end

% function B = vgradmat2(self,N,Ndersp,c)
% nfens= size(Ndersp,1);
% dim=2;
% B = zeros(dim*dim,nfens*dim);
% for i= 1:dim
% B(dim*(i-1)+1:dim*i,i:dim:nfens*dim-dim+i)=Ndersp';
% end
% end

% function B = vgradmat3(self,N,Ndersp,c)
% nfens= size(Ndersp,1);
% dim=3;
% B = zeros(dim*dim,nfens*dim);
% for i= 1:dim
% B(dim*(i-1)+1:dim*i,i:dim:nfens*dim-dim+i)=Ndersp';
% end
% end
