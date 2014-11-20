function  renderLaTeX2E(fontSize)
    % Render the latex string from the clipboard.
    %
    % function  renderLaTeX2E(fontSize)
    %
    % fontSize= font size (24 if not supplied)

    c=clipboardpaste;
    % default parameter
    if(~exist('fontSize'))
        fontSize=24;
    end

    if (~strcmp(c.primaryType,'text'))
        c.data ='\mbox{Clipboard cannot be processed}';
    end

    % Create a figure for rendering the equation.
    hFigure = figure('handlevisibility','off', 'NumberTitle','off',...
        'Units','Characters',...
        'integerhandle','off', ...
        'visible','off', ...
        'paperpositionmode', 'auto', ...
        'color','w');
    set(hFigure,'menubar','none');
    hAxes = axes('position',[0 0 1 1], ...
        'parent',hFigure, ...
        'xtick',[],'ytick',[], ...
        'xlim',[0 1],'ylim',[0 1], ...
        'visible','off');
    hText = text('parent',hAxes,'position',[.5 .5], ...
        'horizontalalignment','center','fontSize',fontSize, ...
        'interpreter','latex');

    % Render and snap!
    [lastMsg,lastId] = lastwarn('');
    c.data = regexprep(c.data, '\n', ' ');
    set(hText,'string', ['$$' c.data '$$']);

    if ispc
        % The font metrics are not set properly unless the figure is visible.
        set(hFigure,'Visible','on');
        set(hFigure,'menubar','none');
    end

    % We adapt the figure size to the equation size, in order to crop it
    % properly when printing. The following lines allow to respect the font
    % size.

    textDimension   =   get(hText,'Extent');
    figureDimension =   get(hFigure,'Position');
    %     p =get(get(hFigure,'CurrentAxes'),'Position');
    % get(hFigure,'Units')
    % get(get(hFigure,'CurrentAxes'),'Units')

    newWidth=textDimension(3)*figureDimension(3) +0.0002;
    newHeight=textDimension(4)*figureDimension(4)+0.0002;

    set(hFigure,'Position',[1 figureDimension(2) newWidth newHeight]);
    %     set(hFigure,'Position',[figureDimension(1) figureDimension(2) newWidth newHeight]);
    set(hFigure,'WindowStyle','modal','ToolBar','none');
    % Draw the figureWindowStyle is modal
    drawnow;
    %   set(get(hFigure,'CurrentAxes'),'Position',[figureDimension(1) figureDimension(2) newWidth newHeight]);

    texWarning = lastwarn;
    lastwarn(lastMsg,lastId)

    set(hFigure,'renderer','Zbuffer');
    %     set(get(gca,'Title'),'fontsize',24)
%     figure(hFigure);
    %     hgexport(hFigure,'-clipboard');
    print(hFigure, '-dbitmap')
    pause(0.001);
    close (hFigure);
end
