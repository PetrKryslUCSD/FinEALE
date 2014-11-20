classdef femm_base
    % Base class for all finite element models
    % from which more specialized (heat diffusion, elasticity,...) finite
    % element model classes are derived.
    %
    %
    %
    
    properties
        material  = [];% material object
        fes = [];% finite element set object
        integration_rule = [];% integration rule object
        Rm=[];% local material orientation matrix
        fen2fe_map = [];% map from finite element nodes to connected finite elements
    end
    
    properties (Access=protected)
        % Associated geometry field. Default is there is none, so any
        % operation that requires a geometry field needs to be supplied it.
        % There may be operations that could benefit from pre-computations
        % that involve a geometry field. If so, associating the geometry
        % field gives the FEMM a chance to save on repeated computations.
        assoc_geom = []; % associated geometry field
    end
    
    methods % constructor
        
        function self = femm_base (Parameters)
            % Constructor.
            % Parameters:
            %  mandatory:
            %      fes = finite element set.  The type of the FE set will be dependent upon
            %        the operations required. For instance, for interior (volume) integrals
            %        such as body load or the stiffness hexahedral H8 may be used whereas
            %        for boundary  (surface) integrals quadrilateral Q4 would be needed.
            %  optional:
            %      mater =material,
            %      integration_rule= integration rule object
            %      Rm= material orientation matrix, or a handle to function to compute the
            %         material orientation matrix, or a string denoting the type 
            %         of material orientation matrix to be used ('geniso').
            %            In the columns of the material orientation matrix are the basis vectors
            %         expressed in the global Cartesian basis. If the orientation matrix
            %         is not supplied, Rm==identity is assumed. The orientation matrix
            %         may be passed in as empty ([]), and the model code can take advantage
            %         of that by assuming that Rm being empty means the local-to-global
            %         transformation is the identity and avoiding thereby the
            %         multiplication.
            %            The function to compute the orientation matrix should have
            %         the signature
            %               function Rm = SampleRm(XYZ, ts, label)
            %         The orientation matrix can be computed based on any of the three
            %         arguments.
            %         XYZ= global Cartesian location of the point at which the orientation
            %                 matrix is desired,
            %         ts= the Jacobian matrix with the tangents to the parametric coordinates
            %                 as columns,
            %         label= the label of the finite element in which the point XYZ resides
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self.material = [];
            self.fes = [];
            self.integration_rule = [];
            self.Rm=[];
            if isfield(Parameters,'material')
                self.material = Parameters.material;
            end
            if isfield(Parameters,'fes')
                self.fes = Parameters.fes;
            end
            if isfield(Parameters,'integration_rule')
                self.integration_rule = Parameters.integration_rule;
            end
            Rm=[];
            if isfield(Parameters,'Rm')
                Rm= Parameters.Rm;
                if (~isempty(Rm)) && strcmp(class(Rm),'double')
                    if norm(Rm'*Rm-eye(size(Rm, 2)))>1e-9
                        error('Non-orthogonal local basis matrix!');
                    end
                elseif strcmp(class(Rm),'function_handle')
                elseif strcmp(class(Rm),'char')
                    switch Rm
                    case 'geniso'
                        Rm=@femm_base.geniso_Rm;
                    otherwise
                        error(['Cannot handle type of Rm: ' Rm])
                    end
                elseif ((~isempty(Rm)))
                    error(['Cannot handle class of Rm: ' class(Rm)])
                end
            end
            self.Rm=Rm;
        end
        
    end
    
    methods
        
        function [npts Ns gradNparams w pc] = integration_data (self)
            % Calculate the data needed for  numerical quadrature.
            %
            % function [npts Ns gradNparams w pc] = integration_data (self)
            %
            integration_rule = self.integration_rule;
            pc = self.integration_rule.param_coords;
            w  =  self.integration_rule.weights ;
            npts = self.integration_rule.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
            end
        end
        
        function Boolean = is_material_orientation_constant(self)
        % Does the material orientation matrix correspond to homogeneous material?
          %
            % Does the material orientation matrix correspond to an homogeneous
            % material  (i.e. is it the same at each location)? If that is not the case
            % we have to keep calling material_orientation() to get the appropriate
            % orientation matrix at each quadrature point.
            %
            % function Boolean = is_material_orientation_constant(self)
            %
            Boolean =false; % assume it is not: Rm would be a function handle
            if isempty(self.Rm)
                Boolean=true;% this implies identity  for the transformation
            else
                if strcmp(class(self.Rm),'double')% the transformation matrix is constant
                    Boolean= true; % constant matrix
                end
            end
        end
        
        function fen2fe_map =node_to_element_map(self)
            % Create a map from nodes to elements.
            %
            % function fen2fe_map =node_to_element_map(self)
            %
            if (isempty(self.fen2fe_map))
                self.fen2fe_map=fenode_to_fe_map (struct ('fes',self.fes));
            end
            fen2fe_map =self.fen2fe_map;
        end
        
        function result = integrate_field_function (self, geom, a_field, fh, varargin)
            % Integrate a nodal-field function over the discrete manifold.
            %
            % function result = integrate_function (self, geom, a_field, fh, varargin)
            %
            % Integrate given function of a given field over the finite elements
            % covering a manifold.
            %
            % Arguments:
            % self = femmlock, geom = geometry field, a_field = an arbitrary field,
            % fh   = function handle, function f(Location,Field_value_at_location)
            % varargin= optional: manifold dimension to be supplied to Jacobian_mdim().
            %
            fes =self.fes;
            [npts Ns gradNparams w pc] = integration_data (self);
            if nargin >=5
                m =varargin{1};% Either the manifold dimension was supplied
            else
                m= fes.dim;% ...Or it is implied
            end
            result = [];
            % Now loop over all fes in the block
            conns = fes.conn; % connectivity
            xs = geom.values;
            for i=1:size(conns,1)
                conn =conns(i,:);
                x = xs(conn,:);
                V = gather_values(a_field, conn);
                % Loop over all integration points
                for j=1:npts
                    J = Jacobian_matrix (fes, gradNparams{j}, x);
                    Jac =Jacobian_mdim(fes, conn, Ns{j}, J, x, m);
                    if (isempty(result))
                        result =fh(Ns{j}'*x,Ns{j}'*V)*Jac*w(j);
                    else
                        result = result + fh(Ns{j}'*x,Ns{j}'*V)*Jac*w(j);
                    end
                end
            end
        end
        
        function result = integrate_function (self, geom, fh, varargin)
            % Integrate a function over the discrete manifold.
            %
            % function result = integrate_function (self, geom, fh, varargin)
            %
            % Integrate some scalar function over the geometric cells. When the scalar
            % function returns just +1 [measure(femm,geom,inline('1'))], the result
            % measures the volume (number of points, length, area, 3-D volume,
            % according to the manifold dimension). When the function returns the mass
            % density, the method measures the mass, when the function returns the
            % x-coordinate equal measure the static moment with respect to the y- axis,
            % and so on.
            %
            % Arguments: self = femmlock, geom = geometry field, fh   = function handle
            % varargin= optional: manifold dimension to be supplied to Jacobian_mdim().
            %
            % Example:
            %     V=measure(femm,geom,@(x)(1))% Volume
            %     S=measure(femm,geom,@(x)(x))% Static moments
            %     CG=S/V% Center of gravity
            %     % Now compute the moments of inertia
            %     I=measure(femm,geom,@(x)(norm(x-CG)^2*eye(3)-(x-CG)'*(x-CG)))
            %     mass=V*rhos;
            %     Inertia=I*rhos;
            fes =self.fes;
            if nargin >=4
                m =varargin{1};% Either the manifold dimension was supplied
            else
                m= fes.dim;% ...Or it is implied
            end
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            [npts Ns gradNparams w pc] = integration_data (self);
            result = 0;% Initialize the result
            % Now loop over all fes in the block
            conns = fes.conn; % connectivity
            xs =geom.values;
            for i=1:size(conns,1)
                conn =conns(i,:);
                x = xs(conn,:);
                % Loop over all integration points
                for j=1:npts
                    J = Jacobian_matrix (fes, gradNparams{j}, x);
                    Jac =Jacobian_mdim(fes, conn, Ns{j}, J, x, m);
                    result = result + fh(Ns{j}'*x)*Jac*w(j); %'
                end
            end
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
            % Default is there is none, so any
            % operation that requires a geometry field needs to be supplied it.
            % There may be operations that could benefit from pre-computations
            % that involve a geometry field. If so, associating the geometry
            % field gives the FEMM a chance to save on repeated computations.
            self.assoc_geom = geom;
        end
        
        function draw(self, gv, context)
            %  Produce graphic representation for all FEs.
            %
            % function draw(self, gv, context)
            %
            %
            % Input arguments
            % self = self
            % gv = graphic viewer
            % context = struct
            % with typically mandatory fields
            %    x=reference geometry field
            %    u=displacement field
            % and with optional fields
            %    facecolor = color of the faces, solid
            %    colorfield =field with vertex colors
            %       only one of the facecolor and colorfield may be supplied
            %    shrink = shrink factor
            % and with any others that a particular implementation of a fe might
            % require or recognize.
            %
            % Call as: draw(hotfemm, gv, struct ('x',geom, 'u',100*geom,...
            %             'colorfield',colorfield, 'shrink',1));
            %
            draw(self.fes, gv, context);
        end
        
        function draw_isosurface(self, gv, context)
            %  Draw graphic representation for iso-surfaces.
            %
            % function draw_isosurface(self, gv, context)
            %
            %
            % See the function draw_mesh()
            %
            draw_isosurface(self.fes, gv, context);
        end
        
        
    end
    
    methods(Static)
    
        function Rm = geniso_Rm(XYZ,tangents,fe_label)
        % Material orientation matrix for isotropic materials.
        %
        % function Rm = geniso_Rm(XYZ,tangents,fe_label)
        %
        % XYZ = location at which the material directions are needed
        % tangents = tangent vectors to parametric coordinates in columns
        % fe_label= label of the finite element
        %
        % The basic assumption here is that the material is isotropic, and
        % therefore the choice of the material directions does not really matter as
        % long as they correspond to the dimensionality of the element. For 
        % instance a one-dimensional element (L2 as an example) may be embedded 
        % in a three-dimensional space.
        %
        % This function assumes that it is being called for
        % an ntan-dimensional manifold element, which is embedded in a
        % sdim-dimensional Euclidean space. If ntan == sdim,
        % the material directions matrix is the identity; otherwise the local
        % material directions are aligned with the linear subspace defined by the
        % tangent vectors.
        % 
        % Warning: this *cannot* be reliably used to produce consistent stresses
        % because each quadrature point gets a local coordinate system which
        % depends on the orientation of the element.
        %
            [sdim, ntan] = size(tangents);
            if sdim==ntan
                Rm=eye(sdim);
            else
                e1=tangents(:,1)/norm(tangents(:,1));
                switch ntan
                    case 1
                        Rm = [e1];
                    case 2
                        n =skewmat(e1)*tangents(:,2)/norm(tangents(:,2));
                        e2=skewmat(n)*e1;
                        e2=e2/norm(e2);
                        Rm = [e1,e2];
                    otherwise
                        error('Got an incorrect size of tangents');
                end
            end
        end

    end

end
