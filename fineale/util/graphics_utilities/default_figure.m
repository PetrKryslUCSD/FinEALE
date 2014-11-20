
function h=default_figure(varargin)
    % Make default figure. 
    %
    %     h=default_figure(varargin)
    %
    % varargin= passed along to figure()
    h=figure(varargin{:});
    set_graphics_defaults(h);