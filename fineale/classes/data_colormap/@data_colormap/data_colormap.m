classdef data_colormap
    % Class of the map from data to colors.
    %
    %
    
    properties
        rmin = [];% minimum of the value range
        rmax = [];% maximum of the value range
        colormap = [];% colormap
    end
    
    methods
        
        function self = data_colormap (Parameters)
            % Constructor.
            % Parameters:
            %    range =array with range_min, and range_max,
            %    colormap =color map, for instance jet, hsv, cadcolors, and so on
            %
            % Example:
            % dcm=data_colormap(struct('range',[min(Uz),max(Uz)],'colormap',jet(16)));
            if nargin < 1
                return
            end
            self.rmin = min(Parameters.range);
            self.rmax = max(Parameters.range);
            if (self.rmax==self.rmin)
                self.rmax=self.rmax+eps(self.rmax);
            end
            
            self.colormap = hot;
            if (isfield(Parameters,'colormap'))
                self.colormap = Parameters.colormap;
            end
        end
        
        function color = map_data (self, v)
            % Map data values to color.
            %
            % function color = map_data (self, v)
            %
            % v= data value, array of numbers to be mapped to colors
            %
            % Call as:
            % c=field(struct ('name', ['c'], 'data', map_data(dcm, Uz)));
            % The field c ('colorfield') may then be passed to graphic drawing primitives.
            %
            n = length(v);
            color = zeros(n,3);
            if (self.rmax > self.rmin)
                Nc=size(self.colormap,1);
                k = ((v-self.rmin)/(self.rmax-self.rmin))*(Nc-2)+1;
                ix=find((k < 1)); 
                if (~isempty(ix)),k(ix)=1;end
                ix=find((k>Nc-1)); 
                if (~isempty(ix)),k(ix)=size(self.colormap,1);end
                color(:,:) = self.colormap(round(k),:);
            end
            return;
        end
        
        function value = range (self)
            % Get the range of the map (minimum, maximum).
            %
            % function value = range (self)
            %
        value =[self.rmin,self.rmax];
        end
        
    end
    
end
