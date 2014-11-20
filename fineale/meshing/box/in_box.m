% Is the given location inside a box?
%
% function b= in_box (box,x)
%
% box = bounding box
%     for 1-D box=[minx,maxx], or
%     for 2-D box=[minx,maxx,miny,maxy], or
%     for 3-D box=[minx,maxx,miny,maxy,minz,maxz] 
% x= location
%
% See also: bounding_box, boxes_overlap, inflate_box,     update_box
%
function b= in_box (box,x)
    dim=length (box)/2;
    b=inrange (box(1:2),x(1));
    for i=2:dim
        b=(b & inrange (box(2*i-1:2*i),x(i)));
    end
end
function b= inrange (range,x)
    b=((x>=range(1)) & (x<=range(2)));
end