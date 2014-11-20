classdef fenode_to_fe_map
    % Constructor of the map from finite element nodes to the finite elements.
    %
    
    properties
        % Map from finite element nodes to the finite elements connecting them.
        % For each  node referenced in the connectivity of 
        % the finite element set on input, the numbers of the individual 
        % finite elements that reference that node is stored in an array in 
        % the so array map{}.
        %         Example: fes.conn= [7,6,5;
        %                             4,1,3;
        %                             3,7,5];
        %             The map reads
        %                     map{1} = [2];
        %                     map{2} = [];%  note that node number 2 is not referenced by the connectivity
        %                     map{3} = [2,3];
        %                     map{4} = [2];
        %                     map{5} = [1,3];
        %                     map{6} = [1];
        %                     map{7} = [1,3];
        % The individual elements from the connectivity that reference 
        % node number 5 are 1 and 3, so that fes.conn(map{5},:) lists all the 
        % nodes that are connected to node 5 (including node 5 itself). 
        map = {};
    end
    
    methods
        
        function self = fenode_to_fe_map (Parameters)
            % Constructor.
            % Parameters:
            %    fes =set of finite elements (descendent of fe_set)
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            self.map = {};
            conns=Parameters.fes.conn;
            n=  max(max(conns));
            [self.map{1:n}] = deal([]);
            cs=size(conns,1);
            for i=1:size(conns,2)
                for j=1:size(conns,1)
                    ni=conns(j,i);
                    self.map{ni}= [self.map{ni}, j];
                end
            end
        end
        
    end
    
end



