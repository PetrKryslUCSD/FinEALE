%Display a status/progress bar.
%
%Display a status/progress bar and inform about the elapsed
%as well as the remaining time (linear estimation).
%
%Synopsis:
%
%  f=progressbar
%     Get all status/progress bar handles.
%
%  f=progressbar(title)
%     Create a new status/progress bar. If title is an empty
%     string, the default 'Progress ...' will be used.
%
%  f=progressbar(title,f)
%     Reset an existing status/progress bar or create a new
%     if the handle became invalid.
%
%  f=progressbar(done,f)
%     For 0 < done < 1, update the progress bar and the elap-
%     sed time. Estimate the remaining time until completion.
%     On user abort, return an empty handle.
%
%  delete(progressbar)
%     Remove all status/progress bars.
%
%  drawnow
%     Refresh all GUI windows.
%
%Example:
%
%     f=progressbar('Wait some seconds ...');
%     for p=0:0.01:1
%        pause(0.2);
%        if isempty(progressbar(p,f))
%           break;
%        end;
%     end;
%     delete(f);
%
%
%        Leutenegger Marcel © September 2004
%        (C) 2006, Petr Krysl
%
function f=progressbar(p,f)
    if nargin < nargout           % get handles
        o='ShowHiddenHandles';
        t=get(0,o);
        set(0,o,'on');
        f=findobj(get(0,'Children'),'flat','Tag','progressbar');
        set(0,o,t);
        return;
    end
    if nargin & ischar(p)
        if nargin == 2 && (~isempty(check(f)))  % reset
            modify(f,'Line','XData',[4 4 4]);
            modify(f,'Rect','Position',[4 54 0.1 22]);
            modify(f,'Done','String','0');
            modify(f,'Time','String','0:00:00');
            modify(f,'Task','String','0:00:00');
        else
            f=create;               % create
        end
        if p
            set(f,'Name',p);
        end
        userdata.t=[cputime cputime 0];
        userdata.Stopped =0;
        set(f,'CloseRequestFcn',@close,'UserData',userdata);
        drawnow;
    elseif nargin == 2 && (~isempty(check(f))) % update
        userdata=get(f,'UserData');
        if isempty(userdata)
            t = [];
        else
        t= userdata.t;
        end
        if any(t < 0)              % confirm
            answer =questdlg({'Stop now? Or resume?',''},...
                    'Action requested','Abort','Stop','Resume','Resume');
            if p >= 1 | strcmp(answer,'Stop')
                userdata.Stopped=~0;
                set(f,'UserData',userdata);
                keyboard;
                userdata.Stopped=0;
            elseif strcmp(answer,'Abort')
                delete(f);
                f= [];
                return
            end
            t=abs(t);
            userdata.t=t;
            set(f,'UserData',userdata);    % continue
        end
        p=min(1,max([0 p]));
        %
        % Refresh display if
        %
        %  1. still computing
        %  2. computation just finished
        %    or
        %     more than a second passed since last refresh
        %    or
        %     more than 0.4% computed since last refresh
        %
        if isempty(t) 
            delete(f);
        end
        if (any(t) & (p >= 1 | cputime-t(2) > 1 | p-t(3) > 0.004))
            userdata=get(f,'UserData');
            userdata.t=[t(1) cputime p];
            set(f,'UserData',userdata);
            t=round(cputime-t(1));
            h=floor(t/60);
            modify(f,'Line','XData',[4 4+504*p 4+504*p]);
            modify(f,'Rect','Position',[4 54 max(0.1,504*p) 22]);
            modify(f,'Done','String',sprintf('%u',floor(p*100+0.5)));
            modify(f,'Time','String',sprintf('%u:%02u:%02u',[floor(h/60);mod(h,60);mod(t,60)]));
            if p > 0.05 | t > 60
                t=ceil(t/p-t);
                h=floor(t/60);
                modify(f,'Task','String',sprintf('%u:%02u:%02u',[floor(h/60);mod(h,60);mod(t,60)]));
            end
            if p == 1
                set(f,'CloseRequestFcn','delete(gcbo);','UserData',[]);
            end
            %             figure(f);%Petr Krysl
            drawnow;
        end
    end
    if ~nargout
        clear;
    end
end

%Check if a given handle is a progress bar.
%
function f=check(f)
    if (~isempty(f)) && ishandle(f(1)) && strcmp(get(f(1),'Tag'),'progressbar')
        f=f(1);
    else
        f=[];
    end
end

%Create the progress bar.
%
function f=create
    s=[512 80];
    t=get(0,'ScreenSize');
    f=figure('DoubleBuffer','on','HandleVisibility','off','MenuBar','none',...
        'Name','Progress ...','IntegerHandle','off','NumberTitle','off',...
        'Resize','off','Position',[floor((t(3:4)-s)/2) s],'Tag','progressbar','ToolBar','none');
    a.Parent=axes('Parent',f,'Position',[0 0 1 1],'Visible','off','XLim',[0 512],'YLim',[0 80]);
    %
    %Horizontal bar
    %
    rectangle('Position',[4 54 504 22],'EdgeColor','white','FaceColor',[0.7 0.7 0.7],a);
    line([4 4 508],[55 76 76],'Color',[0.5 0.5 0.5],a);
    rectangle('Position',[4 54 0.1 22],'EdgeColor','white','FaceColor','green','Tag','Rect',a);
    line([4 4 4],[54 54 77],'Color',[0.2 0.2 0.2],'Tag','Line',a);
    %
    %Description texts
    %
    a.FontWeight='bold';
    a.Units='pixels';
    a.VerticalAlignment='middle';
    text(264,66,1,'%',a);
    text(144,36,'Elapsed time:',a);
    text(144,20,'Remaining:',a);
    text(328,36,'[h:m:s]',a);
    text(328,20,'[h:m:s]',a);
    %
    %Information texts
    %
    a.HorizontalAlignment='right';
    text(264,66,1,'0',a,'Tag','Done');
    text(326,36,'0:00:00',a,'Tag','Time');
    text(326,20,'0:00:00',a,'Tag','Task');
end

%Modify an object property.
%
function modify(f,t,p,v)
    set(findobj(f,'Tag',t),p,v);
end

% close function
function close(src,eventdata)
    userdata =get(gcbo,'UserData');
    if userdata.Stopped
        delete(gcbo);
    else 
        userdata.Stopped=~0;
    end
    userdata.t =-abs(userdata.t);
    set(gcbo,'UserData',userdata);
end