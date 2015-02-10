function box = update_box(box,x)
% Update a box with another location, or create a new box.
% 
% function box = update_box(box,x)
% 
% box = either empty ([]), in which case a box is created using the
% supplied location, or an existing box which is expanded to include the
% supplied location x.   The variable x  can hold multiple points in rows..
% 
% box = bounding box
%     for 1-D box=[minx,maxx], or
%     for 2-D box=[minx,maxx,miny,maxy], or
%     for 3-D box=[minx,maxx,miny,maxy,minz,maxz] 
% See also:  bounding_box, boxes_overlap,  inflate_box,       in_box, 
%
    if isempty(box)
        box=zeros(1,2*size(x,2));
        for i=1:size(x,2)
            box(2*i-1)=Inf;
            box(2*i)=-Inf;
        end
    end
    dim=length(box)/2;
    for j=1:size(x,1)
        for i=1:dim
            box(2*i-1)=min(box(2*i-1),x(j,i));
            box(2*i)=max(box(2*i),x(j,i));
        end
    end
end
