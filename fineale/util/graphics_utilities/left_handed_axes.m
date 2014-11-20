% Set the orientation of a 2-D view for a left-handed coordinate system.
% 
% function left_handed_axes
%     
% Utility function to enable presentation of cable deflections in a
% left-handed coordinate system
function left_handed_axes
    set(gca,'View', [0,-90]);
end
