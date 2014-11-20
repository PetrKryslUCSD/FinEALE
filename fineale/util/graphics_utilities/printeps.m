function printeps(Fig,name, varargin)
% Print the figure to EPS file.
%
% function printeps(Fig,name, varargin)
%
    if (isempty(strfind(name, '.eps')))
        name =[name   '.eps'];
    end
    %     print( Fig,['-r' num2str(resolution)],'-depsc2', name)
    exportfig( Fig,name, 'width',6,'Color', 'rgb', 'fontmode','fixed', varargin{:});
end