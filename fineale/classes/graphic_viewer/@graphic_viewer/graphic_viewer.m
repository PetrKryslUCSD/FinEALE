classdef graphic_viewer
% Graphic viewer class.
%
%

    properties
        axes =[];% Matlab axes handle for drawing
        limits=[];% Limits of the graphics drawn into the graphic viewer
        snap_points =[];% Array of points to which the view can be retargeted (snapped)
    end
    
    properties (Hidden)
        pk_defaults =0;
    end
    
    methods
    
        function self = graphic_viewer (Parameters)
        % Constructor.
        % Parameters: 
        %    axes= Matlab axes handle (optional). If the 'axes' parameter is not
        %        supplied a new figure is created with new axes and those are used by
        %        the graphic viewer.
        %    snap_points = array of points, one per row, to which the view can 
        %        be targeted. (See the rotate() method of the graphic_viewer class.)
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                return
            end
            self.axes =0;
            if (isfield(Parameters,'axes' ))
                self.axes =Parameters.axes;
            else
                h=figure('visible','off');
                self.axes =[];
            end
            self.pk_defaults =0;
            if (isfield(Parameters,'pk_defaults' ))
                self.pk_defaults =Parameters.pk_defaults;
            end
            self.limits=[];
            self.snap_points =[];
            if (isfield(Parameters,'snap_points' ))
                self.snap_points =Parameters.snap_points;
            end
        end
        
        function fig= figure(self)
          try
            child=self.axes;
            while true
              fig=get(child, 'parent' );
              if (strcmpi(get(fig, 'type' ),'Figure')) || (fig==0)
                break;
              end
              child=fig;
            end
          catch
            fig=0;
          end
        end

    end
    
    methods (Hidden)
        function h=gv_default_figure_(varargin)
            h=figure(varargin{:});
        end    
    end
    
end


