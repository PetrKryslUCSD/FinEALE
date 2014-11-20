function abox = inflate_box (box, inflate)
% Inflate the box.
% 
% function abox = inflate_box (box, inflate)
% 
% box=box,
% box = bounding box
%     for 1-D box=[minx,maxx], or
%     for 2-D box=[minx,maxx,miny,maxy], or
%     for 3-D box=[minx,maxx,miny,maxy,minz,maxz] 
% inflate = scalar, amount  by which to increase the box to the left and to the right
%
% See also: bounding_box, boxes_overlap,    in_box,     update_box
%
    abox = box;
    dim=length (box)/2;
    for i=1: dim
        abox(2*i-1)=min(box(2*i-1),box(2*i))- inflate;
        abox(2*i)=max(box(2*i-1),box(2*i))+ inflate;
    end
end
