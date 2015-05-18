function vlist = v_select(v, options)
% Select locations (vertices).
%
% function vlist = v_select(v, options)
% 
% v= array of locations, one location per row
% options = structure with attributes described below
%
% The options structure can be used to specify the location criterion with
% the following options available:
%
% box
% 
% locations are selected because they are inside a box or on its surface.
%
% Example: v_select(v,struct ('box',[1 -1 2 0 4 3])), selects locations
% which are strictly inside the box
%  -1<= x <=1     0<= y <=2     3<= z <=4
%
% cylinder
%
% locations are selected because they are inside a cylinder or on its
% surface. Called as :
%
%   v_select(v,struct ('cylinder',[x, y, R, h]))  for 2D
%   v_select(v,struct ('cylinder',[x, y, z, R, h])) for 3D
%
% the orientation of the cylinder can be changed by supplying the option
% 'orientation'. x,y,z specify the centre of the clyinder, and R and h the
% radius and height.
%
% Example: v_select(v,struct ('cylinder',[0 0 0, 1, 2])), selects locations
% which are strictly inside the cylinder with centre located at (0,0,0),
% having radius 1 and height 2.
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
% the box (or the distance) to make sure some locations which would be on
% the boundary are either excluded or included.
%
% convhull
%
% Select locations inside a covex hull defined by points on its outer
% surface. 
%
% Example:
% 
% To capture all points in a region such as that in Fig. 1
%
%      .         pi/6 radians         .
%       .   _____-----------_____    .
%        .\/_                   _\/ .
%         .     ___--------___     .
%          . ---               ---                                                 
%          |\                    /|                                           
%          | \                  / |                                           
%          |  \                /  |                                           
%          |   \              /   |                                           
%          |    \            /    |                                           
%          |     \          /     /     _                                     
%          \      \        /     /       /\                                     
%           \      \      /     /       /                                      
%            \      \    /     /       /                                       
%             \      \  /     /       /                                        
%              \      \/     /       /                                         
%               \     | (0,0,1)     / h                                           
%                \    |    /       /                                           
%                 \   |   /       /                                            
%                  \  |  /       /                                             
%                   \ | /       /                                              
%                    \|/      \/_                                                
%                      (0,0,0)                                                     
%       Fig. 1: 3D region in which the points                          
%             are to be selected 
%
% The following code would be appropriate
%
%  % create 10 x-y positions along the curved edge
%  h = 1;
%  x,y] = pol2cart(linspace(-pi/3, pi/3, 10)', repmat(h, 10,1));
%
%  % create points at all vertices
%  convhullpoints = [ 0, 0, 0);
%                     x, y, repmat(0, size(x, 1), 1); 
%                     x, y, repmat(1, size(x, 1), 1);  
%                     0, 0, 1); ];
%
%  selection = v_select(v, struct ('convhull', convhullpoints) )
%
%
% Other Options
%
% invert
%
% flag determining whether to invert the selection, defaults to false.
%
% Example: v_select(v,struct ('box', [1 -1 2 0 4 3], 'invert', true))
% selects locations which lie strictly OUTSIDE the box 
%  -1<= x <=1     0<= y <=2     3<= z <=4
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
    elseif isfield(options, 'cylinder')
        inflate = 0;
        if isfield(options,'inflate')
            inflate = (options.inflate);
        end
        dim = size (v,2);
        centre = options.cylinder(1:dim);
        r = options.cylinder(dim+1) + inflate;
        if dim > 2
            h = options.cylinder(dim+2) + inflate;
        else
            h = inf;
        end
        
        if ~isfield(options, 'orientation')
            options.orientation = 'z';
        else
            if ~ischar(options.orientation)
                error('orientation must be a single character, ''x'', ''y'' or ''z''.')
            elseif numel(options.orientation) ~= 1
                error('orientation must be a single character, ''x'', ''y'' or ''z''.')
            end
        end
        if dim < 2
            v = [v, zeros(size (v,1),1)];
            centre = [centre, 0];
        end
        if dim < 3
            v = [v, zeros(size (v,1),1)];
            centre = [centre, 0];
        end
        switch options.orientation
            case 'z'
                vlist = find( incylinder (centre, r+inflate, h+inflate, v(:,1), v(:,2), v(:,3)) );
            case 'x'
                % switch x and z coords
                vlist = find( incylinder (centre, r+inflate, h+inflate, v(:,3), v(:,2), v(:,1)) );
            case 'y'
                % switch y and z coords
                vlist = find( incylinder (centre, r+inflate, h+inflate, v(:,1), v(:,3), v(:,2)) );     
            otherwise
                error('Unrecognised orientation character.')
        end
        nn = size(vlist,1);
    elseif isfield(options, 'convhull')
        % use inhull to find points inside a convex hull defined by the
        % points in options.convhull
        vlist = find(inhull(v, options.convhull));
        nn = size(vlist,1);
    end
    
    
    if isfield (options, 'invert') && options.invert
        % reverse the selection so we select everything except the
        % specified region
        allinds = (1:size (v,1))';
        allinds(vlist) = [];
        vlist = allinds;
        nn = size (vlist,1);
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
function isin = incylinder(centre, r, h, x, y, z)
% determines if points are in a cylinder oriented with its axis along the z
% direction
%
% Syntax
%
% isin = incylinder(centre, r, h, x, y, z)
%
% Input
%
%  centre - 3 element vector containing the x, y and z coodinate of the
%    center of the cylinder
%
%  r - radius of the cylinder
% 
%  h - height, or length of the cylinder along its axis
% 
%  x,y,z - vectors or matrices of the same size containing the x,y and z
%    positions of the coordinates to test
%
% Output
%
%  isin - matrix of the same size as the input x,y and z matrices
%   containing boolean flags, true if the corresponding point lies in the
%   cylinder (or on its surface), or false otherwise.
%

    isbelowz = (z <= (centre(3) + h/2));
    isabovez = (z >= (centre(3) - h/2));
    
    radialdistsquared = (x - centre(1)).^2 + (y - centre(2)).^2;
    
    isinrad = (radialdistsquared <= r^2);
    
    isin = (isbelowz & isabovez & isinrad);

end

