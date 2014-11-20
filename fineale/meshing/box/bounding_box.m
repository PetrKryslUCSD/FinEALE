function box = bounding_box (x)
% Compute the bounding box of a collection of points.
% 
% function box = bounding_box (x)
% 
% x= array of locations, one per row
% Output:
% box=[min(x(:,1)),max(x(:,1))], or
% box=[min(x(:,1)),max(x(:,1)),min(x(:,2)),max(x(:,2))], or
% box=[min(x(:,1)),max(x(:,1)),min(x(:,2)),max(x(:,2)),min(x(:,3)),max(x(:,3))] 
%
% See also: boxes_overlap, inflate_box,     in_box,     update_box
%
    box=zeros(1,2*size(x,2));
    for i=1:size(x,2)
        box(2*i-1)=min(x(:,i));
        box(2*i)=max(x(:,i));
    end
end
