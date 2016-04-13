classdef  fe_set
    % Finite element (FE) set class.
    %
    % The set encloses some part of the computational domain.
    %
    % This class is the base class for FE sets, i.e. every usable FE set has to be
    % derived from it.
    %
    
    properties (GetAccess = public)
        % Connectivity array of the set.
        % The connectivity array lists the nodes that the finite element in 
        % the set connects. The j-th element connects the nodes
        % conn(j,:)
        conn = []; % connectivity array of the set
        % Number of nodes each finite element in this set connects
        % nfens=   number of nodes each  finite element in this finite element set
        %           *should* connect; it is checked
        %           against size(conn,2): if there is a mismatch, an error is reported.
        nfens =0; 
        %         numerical label, supplied for each cell in the set, or
        %           a single number to be applied to all cells
        label = [];%  numerical label of the set, optional
    end
    
    methods % constructor
        
        function retobj = fe_set(Parameters)
            %             Constructor
            % Parameters:
            %     conn=connectivity of the fe_set: numbers of nodes this set connects.
            %  optional
            %     label=numerical label, supplied for each finite element in the set, or
            %           a single number to be applied to all finite elements
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            retobj.nfens = 0;
            if isfield(Parameters,'conn')
                retobj.conn = Parameters.conn;
            end
            if isfield(Parameters,'nfens')
                retobj.nfens = Parameters.nfens;
            end
            if isfield(Parameters,'label')
                retobj.label = Parameters.label; 
                if (length(retobj.label)==1)
                    retobj.label=zeros(size(retobj.conn,1),1)+Parameters.label;
                end
            end
            if (~ isempty(retobj.conn) & (retobj.nfens ~= size(retobj.conn,2)))
                error('Wrong number of nodes!');
            end
        end
        
    end
    
    methods % access methods
        
        function self=set.label(self,val)
% Set the label for the set.
%
% val =either a single value for the entire set, or an array of values one for each connectivity.  
            if (length(val)==size(self.conn,1))
                self.label=reshape(val,size(self.conn,1),1);
            elseif (length(val)==1)
                self.label=val+zeros(size(self.conn,1),1, class(val));
            else
                % No label intended;
                self.label=[];
            end
        end
        
        
        function self=set.conn(self,val)
            %             Get the connectivity of this set.
            
            % Initially nfens is set to zero. It gets then initialized to
            % the correct number. So initially ignore this error.
            if (self.nfens ~= 0) && (self.nfens ~= size(val,2)) 
                error('Wrong number of nodes!');
            end
            self.conn=val;
        end
        
        function val = count(self)
            % Get the number of individual connectivities in the FE set
            %
            % function val = count(self)
            %
            val =size(self.conn,1);
        end
        
    end
    
    methods % modification methods
        
        function self = subset(self,l)
            % Extract a subset of the individual connectivities in the FE set.
            self.conn =self.conn(unique(l),:);
            if (~isempty(self.label))
                self.label=self.label(unique(l),:);
            end
        end
        
        function self = cat(self, other)
            % Concatenate the connectivities of two FE sets
            %
            % function self = cat(self, other)
            % where
            %  self, other = members of the same class, descendent of the fe_set
            self.conn =[self.conn;other.conn];
            if (~isempty(self.label)) || (~isempty(other.label))
                if (isempty(self.label))
                    self.label=zeros(size(self.conn,1),1);
                end
                if (isempty(other.label))
                    other.label=zeros(size(other.conn,1),1);
                end
            end
            self.label=[self.label;other.label];
        end
        
        function self = update_conn(self,NewIDs)
            % Update the connectivity after the IDs of nodes changed.
            %
            % function self = update_conn(self,NewIDs)
            %
            % Inputs:
            % NewIDs= new node IDs. Note that indexes in the conn array point
            % _into_ the  NewIDs array. After the connectivity was updated
            % this is no longer true!
            %
            NewIDs = reshape(NewIDs, 1, []);
            c=self.conn;
            for i=1:size(c,1)
                c(i,:)=NewIDs(c(i,:));
            end
            self.conn=c;
        end
        
    end
    
    methods %(Abstract= true) % pure virtual methods that need to be overridden in specializations
        
        function val = bfun(self, param_coords)
            % Evaluate the basis function matrix.
            % This method is pure virtual, and must be overriden.
            %
            % function val = bfun(self, param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=DIM.
            %
        end
        
        function Nder = bfundpar(self, param_coords)
            % Evaluate the derivatives of the basis function matrix.
            % This method is pure virtual, and must be overriden.
            %
            %
            % function Nder = bfundpar(self, param_coords)
            %
            % Returns an array of NFENS rows, and DIM columns, where
            %    NFENS=number of nodes, and
            %    DIM=number of spatial dimensions.
            % Call as:
            %    Nder = bfundpar(g, pc)
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            % Column j holds derivatives of the basis functions with respect
            % to the parametric coordinate j. Row i holds derivatives of the basis
            % function i with respect to all the parametric coordinates in turn.
            %
        end
        
        function   draw(self, gv, context)
            % Produce a graphic representation of the fe_set
            % This method is pure virtual, and must be overriden.
            %
            %
            % function   draw(self, gv, context)
            %
        end
        
        function   draw_isosurface(self, gv, context)
            % Produce a graphic representation of an isosurface
            % This method is pure virtual, and must be overriden.
            %
            % function   draw_isosurface(self, gv, context)
            %
        end
        
    end
    
    methods % concrete methods
        %
        %         function Ndersp = bfundsp (self, Nder, x)
        %             % Evaluate the spatial derivatives of the basis functions.
        %             %
        %             % function Ndersp = bfundsp (self, Nder, x)
        %             %
        %             %   Call as:    Ndersp = bfundsp (self, Nder, x)
        %             %   where
        %             %      self=FE set whose basis functions derivatives are to be evaluated,
        %             %      Nder=matrix of basis function derivatives wrt parametric
        %             %           coordinates (see bfundpar); size(Nder) = [nbfuns,dim]
        %             %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
        %             %           nbfuns=number of nodes (== number of basis functions) per FE set
        %             %           of the class class(self)
        %             [nbfuns,dim] = size(Nder);
        %             if (size(Nder) ~= size(x))
        %                 error('Wrong dimensions of arguments!');
        %             end
        %             J = x' * Nder;% Compute the Jacobian matrix
        %             Ndersp = Nder / J;% and evaluate the spatial gradients
        %         end
        %
        function J = Jacobian_matrix (self, Nder, x)
            % Evaluate the Jacobian matrix.
            %
            % function J = Jacobian_matrix (self, Nder, x)
            %
            %   where
            %      self=FE set whose basis functions derivatives are to be evaluated,
            %      Nder=matrix of basis function derivatives wrt parametric
            %           coordinates (see bfundpar); size(Nder) = [nbfuns,dim]
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %           nbfuns=number of nodes (== number of basis functions) per FE set
            %           of the class class(self)
            J = x' * Nder;% Compute the Jacobian matrix
        end
        
        function val = in_parametric(self,param_coords)
            % Are given parametric coordinates inside the element parametric domain?
            %
            % function val = in_parametric(self,param_coords)
            %
            %  This function assumes the parametric domain is such that
            %  the basis functions are between 0.0 and 1.0.  If this is not
            %  suitable for a particular geometric cell needs to redefine this
            %  function appropriately.
            %
            % Returns a Boolean: is the point inside, true or false?
            %
            N = bfun(self,param_coords);
            s =intersect(find(N>=0),find(N<=1));
            val = (length(s)==length(N));
        end
        
        function pc = map2parametric(self, x, c, options)
            % Map a spatial location to parametric coordinates.
            %
            % function pc = map2parametric(self, x, c, options)
            %
            % Map a spatial location to parametric coordinates for the geometric cell.
            % self = geometric cell,
            % x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim],
            % c= spatial location (row array),
            % options= Newton's solver control options, such as tolerance and the
            % maximum allowable number of iterations
            %
            % Note: The dimension of the parametric space must match the space in which
            % the spatial coordinates are given.
            %
            % Returns a row array of parametric coordinates if the solution was
            % successful, otherwise NaN are returned.
            %
            ppc = 0*x(1,:)';
            Tolerance = 0.001;
            maxiter =5;
            pc= ppc;
            for i=1:maxiter
                Nder = bfundpar (self, pc);
                J = Jacobian_matrix (self, Nder, x);
                N = bfun(self,pc);
                R = ((N'*x)-c)';
                pc=ppc -J\R;
                if (norm(pc-ppc,inf) <Tolerance)
                    pc=pc';
                    return;
                end
                ppc=pc;
            end
            pc= 0*x(1,:)' +nan;
        end
        
    end
    
end
