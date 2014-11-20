function Boolean = boxes_overlap(box1, box2)
% Do the given boxes overlap?
% 
% function Boolean = boxes_overlap(box1, box2)
% 
% box = bounding box
%     for 1-D box=[minx,maxx], or
%     for 2-D box=[minx,maxx,miny,maxy], or
%     for 3-D box=[minx,maxx,miny,maxy,minz,maxz] 
%
% See also: bounding_box, inflate_box,     in_box,     update_box
%
    Boolean = false;
    dim=length (box1)/2;
    for i=1: dim
        if (box1(2*i-1)>box2(2*i))
            return;
        end
        if (box1(2*i)<box2(2*i-1))
            return;
        end
    end
    Boolean = true;
end
