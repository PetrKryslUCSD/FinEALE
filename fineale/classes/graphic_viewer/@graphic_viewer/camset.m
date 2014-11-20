function camset (self,v)
% Set camera settings for the view. 
%
% function camset (self,v)
%
% Set camera settings for the view. The values to use may be conveniently
% obtained with camget() after the view settings have been adjusted 
% via the user interface.
%
% v = CameraPosition,CameraTarget,CameraUpVector,CameraViewAngle
%     as an array of 10 elements

if (isempty( self.axes ))
        return;
    end
set(self.axes,'CameraPosition',v(1:3));
set(self.axes,'CameraTarget',v(4:6));
set(self.axes,'CameraUpVector',v(7:9));
set(self.axes,'CameraViewAngle',v(10));
end

