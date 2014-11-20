function retobj = interact (self,context)
% Interactively rotate, zoom, dolly, and retarget the view of a 3-D
% plot.
%
% function retobj = interact (self,context)
%
% Arguments:
% context= structure with optional fields:
% snap_points = points in 3-D space to which the camera may be targeted. 
%
% Manipulation of the 3-D view:
% ============================ 
% Double click to restore the original view.
%
% Click with the right mouse button anywhere outside of the displayed
% graphics to get a context menu.  Otherwise the controls are available by
% pressing keys:
%
% hit "z" key over the axes to switch to ZOOM
% hit "r" key over the axes to switch to ROTATION
% hit "o" key over the axes to switch to CAMERA ROLL
% hit "d" key over the axes to switch to DOLLY (sideways translation)
% hit "t" key over the axes to switch to TARGET (setting of the camera target)
% 
% Note that the mode is indicated by the shape of the pointer.
%
% In ROTATION mode:
% press and hold left mouse button to rotate about screen xy axis
% press and hold middle mouse button to rotate about screen z axis
% press '-' to slow down rotation, press '+' to speed up rotation
% In ZOOM mode:
% press and hold left mouse button to zoom in and out
% In DOLLY mode:
% press and hold left mouse button to dolly the camera horizontally and
% vertically
% In TARGET mode:
% Click on the point to which the camera should be targeted. Note: for this
% to work snap_points must be supplied, either to the reset() method or to this method.
% 
    if exist('context') && (~isempty(context)) && isfield(context,'snap_points')
        self.snap_points= context.snap_points;
    end
    drawnow;
    F= find_parent_figure(self.axes);
    set(F,'visible','on'); % we definitely want that figure visible since we want to interact with it
    set(F,'CurrentAxes',self.axes);
    mymanip3d(self.axes,self.snap_points);
    retobj=self;
end

% --------------------------------------------------------------------
% Based upon view3d
% Torsten Vogel 09.04.1999
% tv.volke@bmw.de
% tested under Matlab 5.2
% inspired from rotate3d by The MathWorks, Inc.
% Petr Krysl, March 2011
% -- Introduced setting of camera target
% -- Introduced Variable rotation rate (vdata.Rotation_rate)
% --------------------------------------------------------------------
function mymanip3d(ax,snap_points)
    fig= find_parent_figure(ax);
    viewact(fig,ax,'rot',snap_points);
    return;
    
    
    % ---------------------------------------------- activation ----------
    function viewact(fig,ax,what,snap_points)
        % de-/activates mymanip3d for the given figure
        mymanip3dObj = findobj(allchild(fig),'Tag','mymanip3dObj');
        if strcmp(what,'off')
            if isempty(mymanip3dObj)
                return
            end
            vdata = get(mymanip3dObj,'UserData');
            uirestore(vdata.uistate);
            set(fig,'KeyPressFcn',vdata.oldkeypressfcn)
            delete(mymanip3dObj);
        else
            if isempty(mymanip3dObj)
                mymanip3dObj = makemymanip3dObj(fig,ax,snap_points);
            end
            vdata = get(mymanip3dObj,'UserData');
            vdata.what = what;
            vdata.Rotation_rate=0.1;
            Set_pointer(vdata.what);
            vdata.snap_points=snap_points;
            set(mymanip3dObj,'UserData',vdata);
        end
        hcmenu = uicontextmenu('parent',fig);
        item1 = uimenu(hcmenu, 'Label', 'Zoom mode (z)', 'Callback', @(e,o)mymanip3dContextFcn('z'));
        item1 = uimenu(hcmenu, 'Label', 'Rotate mode (r)', 'Callback', @(e,o)mymanip3dContextFcn('r'));
        item1 = uimenu(hcmenu, 'Label', 'Camera roll mode (o)', 'Callback', @(e,o)mymanip3dContextFcn('o'));
        item1 = uimenu(hcmenu, 'Label', 'Dolly mode (d)', 'Callback', @(e,o)mymanip3dContextFcn('d'));
        item1 = uimenu(hcmenu, 'Label', 'Target mode (t)', 'Callback', @(e,o)mymanip3dContextFcn('t'));
        item1 = uimenu(hcmenu, 'Label', 'Slowdown rotation (-)', 'Callback', @(e,o)mymanip3dContextFcn('-'));
        item1 = uimenu(hcmenu, 'Label', 'Speed up rotation (+)', 'Callback', @(e,o)mymanip3dContextFcn('+'));
        set(fig,'uicontextmenu',hcmenu);
    end
    
    function mymanip3dContextFcn(currchar)
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            if isempty(mymanip3dObj)
                return
            end
            vdata = get(mymanip3dObj,'UserData');
            if strcmp(currchar,'r')
                vdata.what = 'rot';
            elseif strcmp(currchar,'z')
                vdata.what = 'zoom';
            elseif strcmp(currchar,'o')
                vdata.what = 'roll';
            elseif strcmp(currchar,'d')
                vdata.what = 'dolly';
            elseif strcmp(currchar,'t')
                vdata.what = 'target';
                set(gcf,'Pointer','custom','pointershapecdata',pointershapes('target'));
            elseif strcmp(currchar,'+')
                vdata.Rotation_rate= vdata.Rotation_rate*2;
                vdata.Rotation_rate= min([vdata.Rotation_rate,1]);
            elseif strcmp(currchar,'-')
                vdata.Rotation_rate= vdata.Rotation_rate/2;
                vdata.Rotation_rate= max([vdata.Rotation_rate,1e-4]);
            end
            Set_pointer(vdata.what);
            set(mymanip3dObj,'UserData',vdata)
        end
    
    % ---------------------------------------------- make mymanip3dObj ------
    function mymanip3dObj = makemymanip3dObj(fig,ax,snap_points)
        % ---------------------------------------------- mymanip3dDownFcn -------
        function mymanip3dDownFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            mouseclick = get(gcf,'SelectionType');
            if isempty(mymanip3dObj)
                return
            end
            vdata = get(mymanip3dObj,'UserData');
            vdata.oldunits = get(gcf,'Units');
            set(gcf,'Units','pixels');
            vdata.old_pt = get(0,'PointerLocation');
            %  ----------------- store or restore previous view
            ViewData = get(get(gca,'zlabel'),'UserData');
            if isempty(ViewData)
                ViewData = manageViewData('get_from_axes');
                set(get(gca,'zlabel'),'UserData',ViewData)
            end
            if strcmp(mouseclick,'open')
                manageViewData('set_axes',ViewData);
                set(gcf,'Units',vdata.oldunits)
                return
            end
            %  ----------------- display text box
            fig_color = get(gcf,'Color');
            c = sum([.3 .6 .1].*fig_color);
            set(vdata.textbox,'BackgroundColor',fig_color);
            if(c > .5)
                set(vdata.textbox,'Color',[0 0 0]);
            else
                set(vdata.textbox,'Color',[1 1 1]);
            end
            %  ----------------- what to do?
            if strcmp(vdata.what,'rot')
                if strcmp(mouseclick,'normal')
                    set(vdata.textbox,'string','Screen XY Rotation');
                    set(gcf,'WindowButtonMotionFcn',@(s,e)mymanip3dxyFcn);
                    %set(gcf,'Pointer','custom','pointershapecdata',pointershapes('rot'));
                elseif strcmp(mouseclick,'alt')
                    set(vdata.textbox,'string','Screen Z Rotation');
                    set(gcf,'WindowButtonMotionFcn',@(s,e)mymanip3dzFcn);
                    %set(gcf,'Pointer','custom','pointershapecdata',pointershapes('rot'));
                end
            elseif strcmp(vdata.what,'roll')
                set(vdata.textbox,'string','Screen Z Rotation');
                    set(gcf,'WindowButtonMotionFcn',@(s,e)mymanip3dzFcn);
                    %set(gcf,'Pointer','custom','pointershapecdata',pointershapes('rot'));
            elseif strcmp(vdata.what,'dolly')
                set(vdata.textbox,'string','Dolly');
                set(gcf,'WindowButtonMotionFcn',@(s,e)mymanip3ddollyFcn);
                %set(gcf,'Pointer','custom','pointershapecdata',pointershapes('dolly'));
            elseif strcmp(vdata.what,'target')
                p=get(gca,'currentpoint');
                cp=get(gca, 'CameraPosition');
                if (~ isempty(vdata.snap_points))
                    d=diff(p,1);
                    d=d'/norm(d);
                    Distance=zeros(size(vdata.snap_points,1),1);
                    for j=1:size(vdata.snap_points,1)
                        at=vdata.snap_points(j,:)-p(1,:);
                        Distance(j)=norm(at-(at*d)*d');
                    end
                    [C,I] = min(Distance);
                    p=vdata.snap_points(I,:);
                    set(gca,'cameratarget', p);
                end
                vdata.what='rot';
                set(vdata.textbox,'string','Target selected');
                set(gcf,'WindowButtonMotionFcn',[]);
                Set_pointer(vdata.what);
            else
                if strcmp(mouseclick,'normal')
                    set(vdata.textbox,'string','Zoom');
                    set(gcf,'WindowButtonMotionFcn',@(s,e)mymanip3dzoomFcn);
                    %set(gcf,'Pointer','custom','pointershapecdata',pointershapes('zoom'));
                end
            end
            set(mymanip3dObj,'UserData',vdata)
            %             set(vdata.textbox,'visi','on') 
        end
        
        % ---------------------------------------------- mymanip3dUpFcn ---------
        function mymanip3dUpFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            if isempty(mymanip3dObj)
                return
            end
            vdata = get(mymanip3dObj,'UserData');
            if (~isempty(vdata.oldunits))
                set(gcf,'WindowButtonMotionFcn',[],'Units',vdata.oldunits);
            end
            set(mymanip3dObj,'visi','off')
        end
        
        % ---------------------------------------------- mymanip3dkeypressFcn ---
        function mymanip3dkeypressFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            if isempty(mymanip3dObj)
                return
            end
            vdata = get(mymanip3dObj,'UserData');
            currchar = lower(get(gcf,'currentchar'));
            if strcmp(currchar,'r')
                vdata.what = 'rot';
            elseif strcmp(currchar,'o')
                vdata.what = 'roll';
            elseif strcmp(currchar,'z')
                vdata.what = 'zoom';
            elseif strcmp(currchar,'d')
                vdata.what = 'dolly';
            elseif strcmp(currchar,'t')
                vdata.what = 'target';
                set(gcf,'Pointer','custom','pointershapecdata',pointershapes('target'));
            elseif strcmp(currchar,'+')
                vdata.Rotation_rate= vdata.Rotation_rate*2;
                vdata.Rotation_rate= min([vdata.Rotation_rate,1]);
            elseif strcmp(currchar,'-')
                vdata.Rotation_rate= vdata.Rotation_rate/2;
                vdata.Rotation_rate= max([vdata.Rotation_rate,1e-4]);
            end
            Set_pointer(vdata.what);
            set(mymanip3dObj,'UserData',vdata)
        end
        
        % ---------------------------------------------- mymanip3dxyFcn ---------
        function mymanip3dxyFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            vdata = get(mymanip3dObj,'UserData');
            if ( ~isfield(vdata,'Rotation_rate'))
                vdata.Rotation_rate= 0.1;
            end
            new_pt = get(0,'PointerLocation');
            old_pt = vdata.old_pt;
            dx = (new_pt(1) - old_pt(1))*vdata.Rotation_rate;
            dy = (new_pt(2) - old_pt(2))*vdata.Rotation_rate;
            direction = [0 0 1];
            coordsys  = 'camera';
            pos  = get(gca,'cameraposition' );
            targ = get(gca,'cameratarget'   );
            dar  = get(gca,'dataaspectratio');
            up   = get(gca,'cameraupvector' );
            [newPos newUp] = camrotate(pos,targ,dar,up,-dx,-dy,coordsys,direction);
            set(gca,'cameraposition', newPos, 'cameraupvector', newUp);
            vdata.old_pt = new_pt;
            set(mymanip3dObj,'UserData',vdata)
        end
        
        % ---------------------------------------------- mymanip3dzFcn ----------
        function mymanip3dzFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            vdata = get(mymanip3dObj,'UserData');
            new_pt = get(0,'PointerLocation');
            old_pt = vdata.old_pt;
            dy = (new_pt(2) - old_pt(2))*.5;
            camroll(gca,-dy)
            vdata.old_pt = new_pt;
            set(mymanip3dObj,'UserData',vdata)
        end
        
        % ---------------------------------------------- mymanip3dzoomFcn -------
        function mymanip3dzoomFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            vdata = get(mymanip3dObj,'UserData');
            new_pt = get(0,'PointerLocation');
            old_pt = vdata.old_pt;
            dy = (new_pt(2) - old_pt(2))/abs(old_pt(2));
            if (1-dy>0)
                camzoom(gca,1-dy);
                %             tgt=get(gca, 'CameraTarget');
                %             XLim=get(gca, 'XLim');
                %             set(gca, 'XLim',tgt(1)+(XLim-tgt(1))/(1-dy));
                %             YLim=get(gca, 'YLim');
                %             set(gca, 'YLim',tgt(2)+(YLim-tgt(2))/(1-dy));
                %             ZLim=get(gca, 'ZLim');
                %             set(gca, 'ZLim',tgt(3)+(ZLim-tgt(3))/(1-dy));
                vdata.old_pt = new_pt;
                set(mymanip3dObj,'UserData',vdata)
            end
        end
        
        % ---------------------------------------------- mymanip3ddollyFcn --------
        function mymanip3ddollyFcn
            mymanip3dObj  = findobj(allchild(gcf),'Tag','mymanip3dObj');
            vdata = get(mymanip3dObj,'UserData');
            new_pt = get(0,'PointerLocation');
            old_pt = vdata.old_pt;
            %             dx = (new_pt(1) - old_pt(1))/old_pt(1)*1.0;
            %             dy = (new_pt(2) - old_pt(2))/old_pt(2)*1.0; 
            dx = (new_pt(1) - old_pt(1));
            dy = (new_pt(2) - old_pt(2));
            %             tgt=get(g a, 'CameraTarget'); 
            camdolly(gca,-dx,-dy, 0, 'movetarget', 'pixels')
            %             ntgt=get(gca, 'CameraTarget');
            %             XLim=get(gca, 'XLim');
            %             set(gca, 'XLim',XLim-(tgt(1)-ntgt(1)));
            %             YLim=get(gca, 'YLim');
            %             set(gca, 'YLim',YLim-(tgt(2)-ntgt(2)));
            %             ZLim=get(gca, 'ZLim');
            %             set(gca, 'ZLim',ZLim-(tgt(3)-ntgt(3)));
            vdata.old_pt = new_pt;
            set(mymanip3dObj,'UserData',vdata)
        end
        
        % save the previous state of the figure window
        vdata.uistate  = uisuspend(fig);
        % the data structure
        vdata.what     = [];
        vdata.old_pt   = [];
        vdata.textbox  = [];
        vdata.oldunits = [];
        vdata.snap_points=snap_points;
        vdata.oldkeypressfcn = get(fig,'KeyPressFcn');
        % mymanip3dObj
        axposition=get(ax, 'position' );
        mymanip3dObj = annotation('textbox',[0 0 0.2 0.05],'Units','Pixels',...
            'Visible','off','tag','mymanip3dObj');
        vdata.textbox  = mymanip3dObj;
        % store current view
        ViewData = manageViewData('get_from_axes');
        set(get(gca,'zlabel'),'UserData',ViewData);
        % functions
        set(fig,'WindowButtonDownFcn',@(src,eventdata)mymanip3dDownFcn);
        set(fig,'WindowButtonUpFcn',@(src,eventdata)mymanip3dUpFcn);
        set(fig,'WindowButtonMotionFcn',[]);
        set(fig,'ButtonDownFcn',[]);
        set(fig,'KeyPressFcn',@(src,eventdata)mymanip3dkeypressFcn);
        
        
        set(mymanip3dObj,'UserData',vdata);
    end
    % ---------------------------------------------- manage ViewData -----
    function ViewData = manageViewData(how,data)
        if nargin == 1 ; data = [];end
        props = {
            'DataAspectRatio'
            'DataAspectRatioMode'
            'CameraPosition'
            'CameraPositionMode'
            'CameraTarget'
            'CameraTargetMode'
            'CameraUpVector'
            'CameraUpVectorMode'
            'CameraViewAngle'
            'CameraViewAngleMode'
            'PlotBoxAspectRatio'
            'PlotBoxAspectRatioMode'
            'Units'
            'Position'
            'View'
            'Projection'
            };
        if strcmp(how,'get_from_axes')
            ViewData = get(gca,props);
        elseif strcmp(how,'get_stored')
            ViewData = get(get(gca,'zlabel'),'UserData');
        elseif strcmp(how,'set_axes')
            set(gca,props,data)
            ViewData = [];
        end
    end
    % -------------------------------------------------------------------------
    function Set_pointer(what)
        if (strcmp(what,'rot')) || (strcmp(what,'zoom')) || (strcmp(what,'dolly')) || (strcmp(what,'target'))  || (strcmp(what,'roll'))
            set(gcf,'Pointer','custom','pointershapecdata',pointershapes(what));
        else
            set(gcf,'Pointer','arrow');
        end
    end
    % -------------------------------------------------------------------------
    % get some pointer shapes
    function shape = pointershapes(arg)
       if strcmp(arg,'zoom')
            % -- zoom
            shape=[ 2   2   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   1   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN  ;
                2   1   2   1   1   1   1   1   2 NaN NaN NaN   2   2   2   2  ;
                2   1   2   1   1   2   1   1   1   2 NaN   2   1   2   1   2  ;
                2   1   2   1   2 NaN   2   1   1   1   2   1   1   2   1   2  ;
                2   2   2   2 NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   1   1   2  ;
                NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   2   2  ];
        elseif strcmp(arg,'dolly')
            shape=[ NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
                NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
                2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
                2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
                NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
                NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
                NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ];
        elseif strcmp(arg,'rot')
            % -- rot
            shape=[ NaN NaN NaN   2   2   2   2   2 NaN   2   2 NaN NaN NaN NaN NaN ;
                NaN NaN NaN   1   1   1   1   1   2   1   1   2 NaN NaN NaN NaN ;
                NaN NaN NaN   2   1   1   1   1   2   1   1   1   2 NaN NaN NaN ;
                NaN NaN   2   1   1   1   1   1   2   2   1   1   1   2 NaN NaN ;
                NaN   2   1   1   1   2   1   1   2 NaN NaN   2   1   1   2 NaN ;
                NaN   2   1   1   2 NaN   2   1   2 NaN NaN   2   1   1   2 NaN ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
                NaN   2   1   1   2 NaN NaN   2   1   2 NaN   2   1   1   2 NaN ;
                NaN   2   1   1   2 NaN NaN   2   1   1   2   1   1   1   2 NaN ;
                NaN NaN   2   1   1   1   2   2   1   1   1   1   1   2 NaN NaN ;
                NaN NaN NaN   2   1   1   1   2   1   1   1   1   2 NaN NaN NaN ;
                NaN NaN NaN NaN   2   1   1   2   1   1   1   1   1 NaN NaN NaN ;
                NaN NaN NaN NaN NaN   2   2 NaN   2   2   2   2   2 NaN NaN NaN ];
        elseif strcmp(arg,'roll')
            % -- rot
            shape=[
   NaN   NaN   NaN   NaN     2     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN
     2     2     2     2     1     1     1     1     1     1     1     1     2     2   NaN   NaN
     2     1     2    1     1     1     1     1     1     1     1     1    1     2   NaN   NaN
     2     1     1     1    1     2     2     2     2     2     2    1     1     1     2   NaN
     2     1     1     1     2     2   NaN   NaN   NaN   NaN     2     2     1     1     1     2
     2     1     1     1     1     2   NaN   2     2    NaN   NaN     2     2     1     1     2
     2     2     2     2     2     2   2     1     1   2   NaN   NaN     2     1     1     2
   NaN   NaN   NaN   NaN   NaN   2     1     1     1     1   2   NaN     2     1     1     2
     2     2     2     2   NaN   2     1     1     1     1   2   NaN     2     1     1     2
     2     1    1     2   NaN   NaN   2     1     1   2   NaN   NaN     2     1     1     2
     2     1     1     2     2   NaN   NaN  2     2    NaN   NaN     2     2     1     1     2
     2     1     1     1     2     2   NaN   NaN   NaN   NaN     2     2     1     1     1     2
     2     2     1     1    1     2     2     2     2     2     2    1     1     1     2   NaN
   NaN     2     2    1     1     1     1     1     1     1     1     1     1     2     2   NaN
   NaN   NaN     2     2     1     1     1     1     1     1     1     1     2   NaN   NaN   NaN
   NaN   NaN   NaN   NaN     2     2     2     2     2     2     2     2   NaN   NaN   NaN   NaN];
              shape (shape==36)=2;
%               shape (shape~=NaN)=1;0
              
        elseif strcmp(arg,'target')
            shape=[ NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
     2     2     2     2     2     2     2     1     1     2     2     2     2     2     2     2
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN
   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     2   NaN   NaN   NaN   NaN   NaN   NaN];
            
            
        end
    end
end
        
function fig= find_parent_figure(ax)
    child=ax;
    while true
        fig=get(child, 'parent' );
        if (strcmpi(get(fig, 'type' ),'Figure')) || (fig==0)
            break;
        end
        child=fig;
    end
end
