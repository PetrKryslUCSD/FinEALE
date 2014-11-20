% Render figure with anti-aliasing.
% 
% function [varargout] = myaa(varargin)
% 
% MYAA Render figure with anti-aliasing.
%   FIG = MYAA(K,AAMETHOD,FIGMODE,RENDERER)
%   The script renders an anti-aliased version of the current figure,
%   making the graphics look a lot better than in a standard matlab figure.
%   This is useful for publishing results on the web or to better see the
%   fine details in a complex and cluttered plot.
%
%   Parameters:
%     K         integer subsampling factor, default K = 4
%               K is limited to < 17.
%     AAMETHOD  subsampling method used, default METHOD = 'standard'
%               Other: 'standard','noshrink','imresize' (fast),'subpixel'
%     FIGMODE   display mode, default = 'figure'
%               Other: 'figure','flashfig','publish', 'replace'
%     RENDERER  default is the renderer of the current figure.
%               Other: 'opengl','painter','zbuffer', ...
%     FIG       The new anti-aliased figure
%
%   Special cases:
%     1) myaa;              % Fast anti-aliasing of the current figure
%                           % (the original figure is not lost)
%     2) myaa('publish')    % Convenient for publishing Matlab code
%                           % (beware, it kills the original figure)
%   Example 1:
%     spharm2;
%     myaa;
%
%   Example 2 (move the mouse over the anti-aliased figure):
%     xpklein;
%     myaa(4,'imresize','flashfig');
%
%   Example 3:
%     Put the following in test.m
%        %% My test publish
%        %  Testing to publish some antialiased images
%        %
%        spharm2;          % Produce some nice graphics
%        myaa('publish');  % Render an anti-aliased version
%     Then run:
%        publish test.m;
%        showdemo test;
%
%   Example 4 (see all the fine details):
%     line(randn(2500,2)',randn(2500,2)','color','black','linewidth',0.01)
%     myaa(8);
%
%   BUGS:
%     Dotted and dashed lines in plots are not rendered correctly. This is
%     probably due to a bug in Matlab and it will hopefully be fixed in a
%     future version.
%     The OpenGL renderer does not always manage to render an image large
%     enough. Try the zbuffer renderer if you have problems or decrease the
%     K factor.
%
%   See also PUBLISH, SNAPNOW
%
%   (See source code for SNAPNOW for undocumented antialiasing options for 
%    use with publish.)
%
%   Version 1, 2008-08-05
%
%   Author: Anders Brun
%           anders@cb.uu.se
%
function [varargout] = myaa(varargin)

%% Flush any graphics & find out about the current DPI...
drawnow;

nowunits = get(gca,'Units');
set(gca,'Units','pixels');
inpixels = get(gca,'Position'); inpixels = inpixels(3);
set(gca,'Units','inches');
ininches = get(gca,'Position'); ininches = ininches(3);
screen_DPI = inpixels / ininches;

%% Determine the best choice of convolver. If IPPL is available, imfilter
% is much faster. Otherwise it does not matter too much.
% if ippl()%Petr Krysl: ippl not found
%     myconv = @imfilter;
%     mypad = 'same';
% else
    myconv = @conv2;
    mypad = 'valid';
% end

%% Set default options and interpret arguments
if isempty(varargin)
    K = 4;
    try
        imfilter(zeros(2,2),zeros(2,2));
        aamethod = 'imresize';
    catch
        aamethod = 'standard';
    end
    figmode = 'figure';
    renderer = get(gcf,'Renderer');
elseif strcmp(varargin{1},'publish')
    K = 4;
    aamethod = 'noshrink';
    figmode = 'publish';
    renderer = get(gcf,'Renderer');
elseif length(varargin) == 1
    K = varargin{1};
    if K > 16
        error('To avoid excessive use of memory, K has been limited to max 16. Change the code to fix this on your own risk.');
    end
    try
        imfilter(zeros(2,2),zeros(2,2));
        aamethod = 'imresize';
    catch
        aamethod = 'standard';
    end
    figmode = 'figure';
    renderer = get(gcf,'Renderer');
elseif length(varargin) == 2
    K = varargin{1};
    aamethod = varargin{2};
    figmode = 'figure';
    renderer = get(gcf,'Renderer');
elseif length(varargin) == 3
    K = varargin{1};
    aamethod = varargin{2};
    figmode = varargin{3};
    if strcmp(figmode,'publish') && ~strcmp(varargin{2},'noshrink')
        printf('\nThe AAMETHOD was not set to ''noshrink'': Fixed.\n\n');
        aamethod = 'noshrink';
    end
    renderer = get(gcf,'Renderer');
elseif length(varargin) == 4
    K = varargin{1};
    aamethod = varargin{2};
    figmode = varargin{3};
    renderer = varargin{4};
else
    error('Wrong syntax, run: help myaa');
end


%% Capture current figure in high resolution
tempfile = 'temp.png';
current_fig = gcf;
current_renderer = get(current_fig,'Renderer');
set(current_fig,'Renderer',renderer);
current_paperpositionmode = get(current_fig,'PaperPositionMode');
current_inverthardcopy = get(current_fig,'InvertHardcopy');
set(current_fig,'PaperPositionMode','auto');
set(current_fig,'InvertHardcopy','off');
print(current_fig,['-r',num2str(screen_DPI*K)], '-dpng', 'temp.png');
set(current_fig,'InvertHardcopy',current_inverthardcopy);
set(current_fig,'PaperPositionMode',current_paperpositionmode);
set(current_fig,'Renderer',current_renderer);
raw_hires = imread(tempfile);
delete(tempfile);


%% Start filtering to remove aliasing
if strcmp(aamethod,'standard')
    % Subsample hires figure image with standard antialiasing using a
    % butterworth filter
    ff = lpfilter(K*3,K*0.9,2);
    kk = ifftshift(real(ifft2(ff)));
    kk = kk./sum(kk(:));
    a1 = max(min(myconv(single(raw_hires(:,:,1))/(256),kk,mypad),1),0);
    a2 = max(min(myconv(single(raw_hires(:,:,2))/(256),kk,mypad),1),0);
    a3 = max(min(myconv(single(raw_hires(:,:,3))/(256),kk,mypad),1),0);
    raw_lowres = double(cat(3,a1(1:K:end,1:K:end),a2(1:K:end,1:K:end),a3(1:K:end,1:K:end)));
elseif strcmp(aamethod,'subpixel')
    if mod(K,3) ~= 0
        error('K must be a multiple of 3 for subpixel antialiasing');
    end
    % Subsample using subpixel shifts ... similar to ClearType
    ff = lpfilter(K*3,K*0.9,2);
    kk = ifftshift(real(ifft2(ff)));
    kk = kk./sum(kk(:));
    a1 = max(min(myconv(single(raw_hires(:,:,1))/(256),kk,mypad),1),0);
    a2 = max(min(myconv(single(raw_hires(:,:,2))/(256),kk,mypad),1),0);
    a3 = max(min(myconv(single(raw_hires(:,:,3))/(256),kk,mypad),1),0);
    raw_lowres = double(cat(3,a1((K+1)/2:K:end,1:K:end-K),a2((K+1)/2:K:end,1+(K/3):K:end-K),a3((K+1)/2:K:end,1+(2*K/3):K:end)));
elseif strcmp(aamethod,'imresize')
    % This is probably the fastest method available at this moment...
    raw_lowres = single(imresize(raw_hires,1/K,'bilinear'))/256;
elseif strcmp(aamethod,'noshrink')
    % This is very useful for creating plots that the user can resize. It
    % is also used for publishing code, since we have little control over
    % the exact scale that Matlab chooses to export the figures in.
    ff = lpfilter(K*3,K*0.9,2);
    kk = ifftshift(real(ifft2(ff)));
    kk = kk./sum(kk(:));
    a1 = max(min(myconv(single(raw_hires(:,:,1))/(256),kk,mypad),1),0);
    a2 = max(min(myconv(single(raw_hires(:,:,2))/(256),kk,mypad),1),0);
    a3 = max(min(myconv(single(raw_hires(:,:,3))/(256),kk,mypad),1),0);
    raw_lowres = double(cat(3,a1(1:end,1:end),a2(1:end,1:end),a3(1:end,1:end)));
end

%% Place the result in the old figure or in a new one...
if strcmp(figmode,'figure');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    fig = figure;
    sz = size(raw_lowres);
    current_units = get(current_fig,'Units');
    set(current_fig,'Units','pixels');
    pos = get(current_fig,'Position');
    set(current_fig,'Units',current_units);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(figmode,'flashfig');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    fig = figure;
    sz = size(raw_lowres);
    current_units = get(current_fig,'Units');
    set(current_fig,'Units','pixels');
    pos = get(current_fig,'Position');
    set(current_fig,'Units',current_units);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
    set(gcf,'WindowButtonMotionFcn','delete(gcf); set(gcf,''WindowButtonMotionFcn'','''')');
elseif strcmp(figmode,'replace');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    set(current_fig,'Units','pixels');
    pos = get(current_fig,'Position');
    close(current_fig);
    fig = figure(current_fig);
    sz = size(raw_lowres);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(figmode,'publish');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    fig = figure;
    current_units = get(current_fig,'Units');
    set(current_fig,'Units','pixels');
    pos = get(current_fig,'Position');
    set(current_fig,'Units',current_units);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    image(raw_lowres);
    set(ax,'Units','normalized');
    set(ax,'Position',[0 0 1 1]);
    axis off;
    close(current_fig);
end


%% Avoid unnecessary console output
if nargout == 1
    varargout(1) = {fig};
end

%% A lowpass filter.
% sz is the size of the filter
% subsmp is the downsampling factor to be used later
% n is the degree of the butterworth filter
function f = lpfilter(sz, subsmp, n)
cutfreq = 0.5 / subsmp;
if mod(sz,2)
    range = (-(sz-1)/2:(sz-1)/2)/(sz-1);
else
    range = (-sz/2:(sz/2-1))/sz;
end
[xx,yy] = ndgrid(range,range);
rr = sqrt(xx.^2+yy.^2);
f = ifftshift(1./(1+(rr./cutfreq).^(2*n)));

