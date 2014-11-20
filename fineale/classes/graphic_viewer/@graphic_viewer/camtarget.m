function camtarget (self,v)
% Set camera  target for the view. 
%
% function camtarget(self,v)
%
% Set CameraTarget settings for the view. 
% 
if (isempty( self.axes ))
    return;
end
set(self.axes,'CameraTarget',v);
end

