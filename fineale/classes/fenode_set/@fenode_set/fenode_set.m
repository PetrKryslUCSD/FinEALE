classdef fenode_set
    % Finite element node set class.
    %
    % function self = fenode_set (varargin)
    %
    % Class representing finite element node set, in one dimension,
    % two dimensions, or three dimensions.
    %
    %
    %
    properties (GetAccess = public)
    % Array of node locations.
    % Array of coordinates, the number of rows corresponds to the number of 
    % nodes in the set and the columns corresponds to the space dimensions.
      % The location of node j is given by x(j,:)
        xyz = [];
    end
    
    methods
        
        function self = fenode_set (Parameters)
        % Constructor.
        % Parameters:
    % xyz= coordinates as an array of dimension [nfens,dim].
    %      nfens= number of nodes, and dim=number of space dimensions 
    %      One row per node and the nodes are indexed from 1 to size(xyz,1).  
    %      So the node numbers are equal to the row numbers.
    % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
        if (nargin ==1)
                self.xyz = Parameters.xyz;
            end
        end
        
    end
    
    methods
        
        function val = dim(self)
        % Get the  dimension of the coordinate that defines the location  of the node.
            val =size(self.xyz,2);
        end
        
        function val= xyz3 (self)
        % Get the  3-D coordinate that define the location  of the node.
        % Even if the nodes  were specified in  lower dimension (1-D, 2-D)
        % this function returns  a 3-D coordinate  by padding with zeros.
            if (size(self.xyz,2)==1)
                val = [self.xyz zeros(size(self.xyz,1),2)];
            elseif (size(self.xyz,2)==2)
                val = [self.xyz zeros(size(self.xyz,1),1)];
            else
                val = self.xyz;
            end
        end
        
        function val = nfens(self)
        % Get the number of finite element nodes in the node set.
            val =count(self);
        end
        
        function val = count(self)
        % Get the number of finite element nodes in the node set.
            val =size(self.xyz,1);
        end
        
        function draw (self, gv, context)
        % Produce a graphic representation of the node set.
        %
        % function draw (self, gv, context)
        %
        %
        % Input arguments
        % self = self
        % gv = graphic viewer
        % context = struct
        % with mandatory fields
        %    x=reference geometry field
        %    u=displacement field
        % with optional fields
        %    color = color of the text
        %
            for i=1:size(self.xyz,1)
                conn = i; 
                x = gather_values(context.x, conn); % coordinates of nodes
                u = gather_values(context.u, conn); % coordinates of nodes
                draw_text (gv, x+u, num2str(i), context);
            end
        end
        
    end
    
end



