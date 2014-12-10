function retobj = clear (self, options)
% "Clear" the associated axes by removing all graphics. 
%
% function retobj = clear (self, options)
%
% options= currently not used
% 
% This method deletes from the Associated axes all graphics objects 
% and resets all axes properties, except Position and Units, 
% to their default values.
    if (isempty(self.axes))||(self.axes==0)
        self.axes=get(gcf,'CurrentAxes');
        if (isempty(self.axes))
            axes;
            self.axes=gca;
        end
        set(self.axes,'Units','normalized')
        set(self.axes,'position',[0.05  0.1  0.8  0.8]);
    end
    try
        %         Parent =get(self.axes, 'parent');
        %         Children =get(Parent,'Children');
        %         for i=1:length(Children)
        %             if (strcmp(get(Children(i), 'type'),'axes'))
        %                  set(Parent,'CurrentAxes',Children(i),'visible',get(Children(i),'visible')); % set the current axes
        %         if (strcmp(get(gcf,'visible'),'on'))
        %             figure(gcf); clf;
        %         else
        %              clf(gcf); self.axes =[];
        %         end
        if  (~isempty(self.axes))&&(strcmp(get(self.axes,'visible'),'on'))
            cla(self.axes,'reset');
        else
            self.axes=get(gcf,'CurrentAxes');
            cla(self.axes); set(self.axes,'visible','off');
        end
%                  break;
%             end
%         end
    catch
        figure(gcf);
        axes;
        self.axes=gca;
    end
    retobj=self;
end

% Author: 	Michael Kleder
% FIGUREC - create a figure window in a non-overlapping (cascading)
%           location
%
% USAGE:
%
% figurec
% figurec(...)
% h=figurec
% h=figurec(...)
%
% FIGUREC acts just like the Matlab FIGURE command, with all arguments
% passed through, except that the new figure is created a little to the
% right and down from the highest numbered figure currently existing, so
% that they won't overlap. If moving the location would push the figure too
% close to the edge of the screen, then the new figure is created in the
% default location as usual. (Subsequent figures will again be cascaded.)
%
% EXAMPLE:
%
% close all;for n=1:20;figurec('color',rand(1,3));plot(0,0);title('Sample');end

function varargout=figurec(varargin)
    f = findobj(0,'type','figure'); % list of existing figures
    ss=get(0,'ScreenSize'); % pixel size of entire screen
    h=figure(varargin{:}); % create figure using pass-through arguments
    hp = get(h,'pos'); % size of new figure when created
    if ~isempty(f)
        f=max(f);
        u=get(f,'units');
        set(f,'units','pixels')
        p=get(max(f),'pos');
        set(f,'units',u)
        % if moving won't push too far, move; else leave in default location
        if p(1)+50+hp(3) <= ss(3) & p(2) >= 5
            u=get(h,'units');
            ss = get(0,'screensize');
            set(h,'units','pixels')
            set(h,'pos',[p(1)+50 p(2)-50 hp(3:4)]);
            set(h,'units',u)
        end
    end
    if nargout > 0
        varargout{1}=h;
    end
end
