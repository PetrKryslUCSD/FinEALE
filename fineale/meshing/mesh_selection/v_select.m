function vlist = v_select(v, options)
% Select locations (vertices).
%
% function vlist = v_select(v, options)
% 
% v= array of locations, one location per row
% options = structure with attributes described below
%
% box
% 
% Select locations using some criterion, for instance and location is selected
% because it is inside a box or on its surface.
% Example: v_select(v,struct ('box',[1 -1 2 0 4 3])), selects
% locations which are strictly inside the box  
%  -1<= x <=1     0<= y <=2     3<= z <=4
% 
% distance
% 
% Example: v_select(v,struct ('distance',0.5, 'from',[1 -1])), selects
% nodes which are less than 0.5 units removed from the point [1 -1].
%
% plane
% 
% Example: v_select(v,struct ('plane',[[1, 0.5, 0.2], -2.3], 'thickness',0.5,')), selects
% nodes which are less than 0.5 units removed from the plane with normal [1, 0.5, 0.2], 
% (the normal is assumed to be of unit length, if it isn't as supplied, 
% it will be normalized internally), at signed distance -2.3 from the
% origin.
%
% nearestto
% 
% Example: v_select(v,struct ('nearestto',[1 -1])), selects
% the node nearest to the point [1 -1].
%
% The option 'inflate' may be used to increase or decrease the extent of
% the box (or the distance) to make sure some locations which would be on the
% boundary are either excluded or included.
% 
    vlist= zeros(1,size(v,1)); nn= 0;
    if isfield(options,'box')
        box=options.box;
        inflate =0;
        if isfield(options,'inflate')
            inflate = (options.inflate);
        end
        dim=length (box)/2;
        for i=1: dim
            abox(2*i-1)=min(box(2*i-1),box(2*i))- inflate;
            abox(2*i)=max(box(2*i-1),box(2*i))+ inflate;
        end
        for i=1:size(v,1)
            if inbox (abox,v(i,:))
                nn =nn +1; vlist(nn) =i;
            end
        end
    elseif isfield(options,'distance')
        from =0*v(1,:);
        if isfield(options,'from')
            from = options.from;
        end
        inflate =0;
        if isfield(options,'inflate')
            inflate = (options.inflate);
        end
        d=options.distance+inflate;
        for i=1:size(v,1)
            if norm (from-v(i,:))<d
                nn =nn +1; vlist(nn) =i;
            end
        end
    elseif isfield(options,'plane')
        n = options.plane(1:3);
        n=n/norm(n);
        inflate =0;
        if isfield(options,'inflate')
            inflate = (options.inflate);
        end
        t=options.thickness+inflate;
        d=options.plane(4);
        for i=1:size(v,1)
            ad=dot(v(i,:),n);
            if abs(d-ad)<t
                nn =nn +1; vlist(nn) =i;
            end
        end
    elseif isfield(options,'nearestto')
        location = options.nearestto;
        distance = ipdm(v, location);
        % Get array index of smallest distance to location from
        % all nodes
        [junk,IX] = min(distance);
        vlist = [IX(1)];
        nn=1;
    end
    if (nn==0)
        vlist = [];
    else
        vlist =vlist(1:nn);
    end
end

function b= inrange (range,x)
    b=((x>=range(1)) && (x<=range(2)));
end
function b= inbox (box,x)
    dim=length (box)/2;
    b=inrange (box(1:2),x(1));
    for i=2:dim
        b=(b && inrange (box(2*i-1:2*i),x(i)));
    end
end

