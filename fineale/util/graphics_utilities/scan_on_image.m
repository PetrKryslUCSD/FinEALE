% Collects points by scanning a curve in a graph image.
%
% function xy = scan_on_image(image_name,x_range,y_range)
%
function xy = scan_on_image(image_name)
    img = imread(image_name);
    image(img,'CDataMapping','scaled');
    axis image ij;
    grid on;
    hold on
    %Pick out the lower left corner of the axes
    x_range =input('x-range?');
    y_range =input('y-range?');
    disp(['Click with left mouse button to pick lower-left corner at ' num2str([x_range(1),y_range(1)])])
    [xi,yi,but] = ginput(1);
    ll=[xi,yi];
    %Pick out the upper right corner of the axes
    disp(['Click with left mouse button to pick upper-right corner at ' num2str([x_range(2),y_range(2)])])
    [xi,yi,but] = ginput(1);
    ur=[xi,yi];
    xy = [];
    % Loop, picking up the points.
    disp('Left mouse button picks points.')
    disp('Right mouse button picks last point.')
    n = 0;
    while true
        [xi,yi,but] = ginput(1);
        %         xi,yi
        if but ~= 1
            break;
        else
            plot(xi,yi,'ro')
            n = n+1;
            xy(:,n) = [xi;yi];
        end
    end
    xy= [x_range(1)+(xy(1,:)'-ll(1))/(ur(1)-ll(1))*diff(x_range) ...
        y_range(1)+ (xy(2,:)'-ll(2))/(ur(2)-ll(2))*diff(y_range)];
end
