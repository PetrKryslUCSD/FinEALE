classdef femm_deformation_linear_nice < femm_deformation_linear
    % Class for the small displacement, small strain deformation
    % model for Nodally-Integrated Continuum Elements (NICE).
    %
    % The approximation is  originally from Dohrmann et al IJNME 47 (2000).
    % The formulation was subsequently developed in Krysl, P. and Zhu, B.
    % Locking-free continuum displacement finite elements with nodal
    % integration, International Journal for Numerical Methods in Engineering,
    % 76,7,1020-1043,2008.
    %
    % See discussion of constructors in <a href="matlab:helpwin 'FAESOR/classes/Contents'">classes/Contents</a>.
    % Options: those recognized by fem_deformation_linear plus
    %
    
    properties (Hidden, SetAccess = private)
        hadjugate = [];
        bfun_gradients = [];
    end
    
    properties
        stabfact=0;
        hBlmat_deviatoric=[];
        hPhi=[];
    end
    
    methods % constructor
        
        function self = femm_deformation_linear_nice (Parameters)
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_deformation_linear(Parameters);
            % We need an appropriate method to compute the adjugate matrix,
            % and we need to override the strain-displacement matrix method
            switch self.fes.dim
                case 0
                    self.hadjugate=[];
                case 1
                    self.hadjugate=[];
                case 2
                    if self.fes.axisymm
                        self.hadjugate=[];
                    else
                        self.hadjugate=@adjugate2;
                        self.hBlmat=@Blmat2_nice;
                    end
                case 3
                    self.hadjugate=@adjugate3;
                    self.hBlmat=@Blmat3_nice;
                    self.hBlmat_deviatoric=@Blmat3_deviatoric_nice;
                    self.hPhi=@Phi3;
                otherwise
                    error (' Not implemented');
            end
            if (isfield(Parameters,'stabfact'))
                self.stabfact=Parameters.stabfact;
            end
        end
        
    end
    
    methods (Hidden)
        
        function p = patch_conn (self,kconns,thisnn)
            % Generate patch connectivity for a given node (thisnn)
            % from the connectivities of the finite elements attached to it.
            %
            %    function p = patch_conn (self,kconns,thisnn)
            %
            p = [setdiff(unique(reshape(kconns,prod(size(kconns)),1)),thisnn)',thisnn];
        end
        
        function bfun_gradients = nodal_bfun_gradients (self, geom)
            % Compute the nodal basis function gradients.
            % Return the cell array of structures with attributes
            %		 bfun_gradients{nix}.Nspd= basis function gradient matrix
            %        bfun_gradients{nix}.Vpatch= nodal patch volume
            %        bfun_gradients{nix}.patchconn= nodal patch connectivity
            %
            
            % If possible, that is if there is an associated geometry, we
            % already have the basis function gradients.  We just return
            % the ones that were precomputed.
            if (~isempty(self.assoc_geom)) && (~isempty(self.bfun_gradients))
                bfun_gradients= self.bfun_gradients;
                return
            end
            fes = self.fes;% grab the finite elements to work on
            % Integration rule: compute the data needed for  numerical quadrature
            [npts Ns Nders w] = integration_data (self);
            % Material orientation
            Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
            if (~Rm_constant)
                Rmh = self.Rm;% handle to a function  to evaluate Rm
            else
                Rm = self.Rm;
            end
            % Material
            mat = self.material;
            % Retrieve data for efficiency
            conns = fes.conn; % connectivity
            labels = fes.label; % connectivity
            xs =geom.values;
            % Get the inverse map from finite element nodes to geometric cells
            fen2fe_map =node_to_element_map(self);
            gmap=fen2fe_map.map;
            % Initialize the nodal gradients, nodal patch, and patch connectivity
            bfun_gradients= cell(length(gmap),1);
            for nix=1:length(gmap)
                bfun_gradients{nix}.Nspd=[];
                bfun_gradients{nix}.Vpatch=0;
                bfun_gradients{nix}.patchconn=[];
            end
            % Now loop over all finite element nodes
            lnmap =0*(1:length(gmap));
            patchconn = [];
            iem =1;
            for nix=1:length(gmap)
                gl=gmap{nix};
                thisnn=nix;
                if ~isempty(gl) % This node has a element patch in this block
                    % establish local numbering of all nodes of the patch @node nix
                    kconns = conns(gl,:); % connectivity
                    patchconn = patch_conn (self,kconns,thisnn);
                    np =length(patchconn);
                    lnmap(patchconn) =1:np;% now store the local numbers
                    c = xs(thisnn, :); % coordinates of central node
                    if (~Rm_constant)% do I need to evaluate the local material orientation?
                        Rm =Rmh(c,[],thisnn);
                    end
                    Nspdavg=zeros(np,fes.dim);% preallocate strain-displacement matrix
                    %                                         Nspdavg=sym(Nspdavg);%!  Symbolic only
                    Vpatch_dim =0; Vpatch =0;
                    for k=1:length(gl)
                        kconn =kconns(k,:);
                        pci=find(kconn == thisnn);% at which node in the element are we?
                        pci =pci(1);% in case there are repeated nodes
                        if (pci<=npts)
                            px = xs(kconn, :); % coordinates of nodes
                            lx=px-ones(length(kconn),1)*c;% centered
                            if (isempty(Rm)),  xl=lx;
                            else               xl=lx*Rm;% local coordinates
                            end
                            J = Jacobian_matrix(fes,Nders{pci},xl);
                            Jac = Jacobian_volume(fes,kconn, Ns{pci}, J, xl);
                            Vpatch=Vpatch+Jac*w(pci);
                            Jac_dim = Jacobian_mdim(fes, kconn, Ns{pci}, J, xl, fes.dim);
                            Nspd = Nders{pci}*self.hadjugate(self,J);
                            Vpatch_dim=Vpatch_dim+Jac_dim*w(pci);
                            Nspdavg(lnmap(kconn),:)=Nspdavg(lnmap(kconn),:)+w(pci)*Nspd;
                        end
                    end
                    if (Vpatch_dim~=0)
                        Nspdavg=Nspdavg/ Vpatch_dim;
                    else
                        Nspdavg=0*Nspdavg;
                        %                         error('Negative nodal volume')
                    end
%                     if (Vpatch<0) %!  Symbolic only
%                         error
%                         Vpatch=0;% If negative volume, ignore it.   However,  at some point we might wish to check for this and raise an error.
%                     end
                    bfun_gradients{nix}.Nspd=Nspdavg;
                    bfun_gradients{nix}.Vpatch=Vpatch;
                    bfun_gradients{nix}.patchconn=patchconn;
                end
                lnmap(patchconn) = 0;
            end
        end
        
        % Finds the adjugate of square matrix 3 x 3 A
        function B = adjugate3(self,A)
            % A = rand(3)
            B =0*A;
            B(1,1) =(A(5)*A(9)-A(8)*A(6));
            B(1,2) =-(A(4)*A(9)-A(7)*A(6));
            B(1,3) =(A(4)*A(8)-A(7)*A(5));
            
            B(2,1) =-(A(2)*A(9)-A(8)*A(3));
            B(2,2) =(A(1)*A(9)-A(7)*A(3));
            B(2,3) =-(A(1)*A(8)-A(7)*A(2));
            
            B(3,1) =(A(2)*A(6)-A(5)*A(3));
            B(3,2) =-(A(1)*A(6)-A(4)*A(3));
            B(3,3) =(A(1)*A(5)-A(4)*A(2));
            % det(A)*inv(A)-B
        end
        
        % Finds the adjugate of square matrix  2 x 2 A
        function B = adjugate2(self,A)
            %     A = rand(2)
            B =0*A;
            B(1,1) =A(2,2);
            B(1,2) =-A(1,2);
            
            B(2,1) =-A(2,1);
            B(2,2) =A(1,1);
            %     det(A)*inv(A)-B
        end
        
        function B = Blmat3_nice(self,Ndersp,c,Rm)
            nfn= size(Ndersp,1);
            B = zeros(6,nfn*3); %initialize
            %                         B =sym(B);%!  Symbolic only
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
        end
        
        function B = Blmat3_deviatoric_nice(self,Ndersp,c,Rm)
            nfn= size(Ndersp,1);
            B = zeros(6,nfn*3); %initialize
            if (isempty(Rm)) % there is no global-to-local transformation
                for i= 1:nfn
                    B(:,3*(i-1)+1:3*i)...
                        = [ (2/3)*Ndersp(i,1) (-1/3)*Ndersp(i,2) (-1/3)*Ndersp(i,3)  ; ...
                        (-1/3)*Ndersp(i,1) (2/3)*Ndersp(i,2) (-1/3)*Ndersp(i,3)  ; ...
                        (-1/3)*Ndersp(i,1) (-1/3)*Ndersp(i,2) (2/3)*Ndersp(i,3) ; ...
                        Ndersp(i,2) Ndersp(i,1) 0  ; ...
                        Ndersp(i,3) 0           Ndersp(i,1) ; ...
                        0           Ndersp(i,3) Ndersp(i,2) ];
                end
            else % global-to-local transformation is requested
                RmT=Rm';
                for i= 1:nfn
                    B(:,3*(i-1)+1:3*i)...
                        = [ (2/3)*Ndersp(i,1) (-1/3)*Ndersp(i,2) (-1/3)*Ndersp(i,3)  ; ...
                        (-1/3)*Ndersp(i,1) (2/3)*Ndersp(i,2) (-1/3)*Ndersp(i,3)  ; ...
                        (-1/3)*Ndersp(i,1) (-1/3)*Ndersp(i,2) (2/3)*Ndersp(i,3) ; ...
                        Ndersp(i,2) Ndersp(i,1) 0  ; ...
                        Ndersp(i,3) 0           Ndersp(i,1) ; ...
                        0           Ndersp(i,3) Ndersp(i,2) ]*RmT;
                end
            end
        end
        
        function B = Blmat2_nice(self,Ndersp,c,Rm)
            nfn= size(Ndersp,1);
            dim =size(Ndersp,2);
            B = zeros(3,nfn*dim);
            if (isempty(Rm)) % there is no global-to-local transformation
                for i= 1:nfn
                    B(:,2*(i-1)+1:2*i)=...
                        [Ndersp(i,1) 0; ...
                        0           Ndersp(i,2); ...
                        Ndersp(i,2) Ndersp(i,1) ];
                end
            else % global-to-local transformation is requested
                RmT=Rm(:,1:2)';
                for i= 1:nfn
                    B(:,2*(i-1)+1:2*i)=...
                        [Ndersp(i,1) 0; ...
                        0           Ndersp(i,2); ...
                        Ndersp(i,2) Ndersp(i,1) ]*RmT;
                end
            end
        end
        
        function Phi =Phi3(self,dim,np,lx,c)
            lx(:,1)=lx(:,1)-c(1);
            lx(:,2)=lx(:,2)-c(2);
            lx(:,3)=lx(:,3)-c(3);
            Phi=zeros(dim*np, 12);
            Phi(1:dim:end-1, 1)=1;
            Phi(1:dim:end-1, 2)=lx(:,1)';
            Phi(1:dim:end-1, 3)=lx(:,2)';
            Phi(1:dim:end-1, 4)=lx(:,3)';
            Phi(2:dim:end, 5)=1;
            Phi(2:dim:end, 6)=lx(:,1)';
            Phi(2:dim:end, 7)=lx(:,2)';
            Phi(2:dim:end, 8)=lx(:,3)';
            Phi(3:dim:end, 9)=1;
            Phi(3:dim:end, 10)=lx(:,1)';
            Phi(3:dim:end, 11)=lx(:,2)';
            Phi(3:dim:end, 12)=lx(:,3)';
        end
        
        
        function self=associate_geometry(self,geom)
            % Associate a geometry field with the FEMM.
            %
            %         function self=associate_geometry(self,geom)
            %
            % geom= Geometry field to associate with self. Pass in []
            %    (empty array) to disassociate any existing geometry from the
            %    self FEMM.
            %
            % By default no geometry field is associated, so any
            % operation that requires a geometry field needs to be supplied it.
            %
            % There may be operations that could benefit from pre-computations
            % that involve a geometry field. If so, associating the geometry
            % field gives the FEMM a chance to save on repeated computations.
            % However, care must be taken after the geometry field was
            % associated: if any of the operations is supplied a different
            % geometry field from the one associated with the FEMM, the
            % result of the operation will be wrong.
            self.bfun_gradients = nodal_bfun_gradients (self, geom);
            self = associate_geometry@femm_base(self, geom);
        end
        
        
    end
    
end
