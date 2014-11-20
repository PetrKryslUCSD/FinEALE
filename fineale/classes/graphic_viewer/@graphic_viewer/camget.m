function v = camget (self)
% Get camera settings for the view.
%
% function v = camget (self)
%
%Output: array of 10 numbers:
% CameraPosition,CameraTarget,CameraUpVector,CameraViewAngle
    v=[];
    if (isempty( self.axes ))
        return;
    end
    v(1:3)=get(self.axes,'CameraPosition');
    v(4:6)=get(self.axes,'CameraTarget');
    v(7:9)=get(self.axes,'CameraUpVector');
    v(10)=get(self.axes,'CameraViewAngle');
end

