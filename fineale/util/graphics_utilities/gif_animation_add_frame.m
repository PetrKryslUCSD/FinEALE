function gif_animation_add_frame(f,frame,Moviefilename,delay,varargin)
% Add a frame to gif animation.
%
% function gif_animation_add_frame(f,frame,Moviefilename,delay,varargin)
%
% f= figure
% frame= frame number; if equal to 1 then the animated gif should be 
%          initialized, otherwise the frame is appended.
% Moviefilename= name of the animated gif
% delay= delay between frames, refer to the argument 'delaytime' of imwrite()
% varargin= these arguments are passed along to  imwrite()
    A = getframe(f);
    [IND, map] = RGB2IND256(A.cdata(:,:,:));
    if frame==1
        imwrite(IND,map,Moviefilename,'gif','WriteMode','overwrite','delaytime',delay/1000,'LoopCount',inf, varargin{:});
    else
        imwrite(IND,map,Moviefilename,'gif','WriteMode','append','delaytime',delay/1000,varargin{:});
    end
end

%RGB2IND256 Naive true color to indexed color converter
function [z,map]=RGB2IND256(rgb)
[r,g,b]=ndgrid(linspace(0,255,8),linspace(0,255,8),linspace(0,255,4));
map=round([r(:) g(:) b(:)])/255;
z=uint8(bitshift(bitand(rgb(:,:,1),224),-5)+bitshift(bitand(rgb(:,:,2),224),-2)+bitand(rgb(:,:,3),192));
end
