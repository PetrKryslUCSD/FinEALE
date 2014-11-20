function felist = fe_select(fens, fes, options)
% Select finite elements.
%
% function felist = fe_select(fens, fes, options)
%
% Select finite elements (FEs) using some criterion, for instance a fe is
% selected because all its nodes are inside a box.
%
% Examples of selection criteria:
%
% box
%
% Select all FEs with all nodes inside the box:
%     fe_select(fens,fes,struct ('box',[1 -1 2 0 4 3]))
% Select all FEs with at least one node inside the box:
%     fe_select(fens,fes,struct ('box',[1 -1 2 0 4 3],'anynode',true))
%
% label
%
% Select all FEs with given label:
%     fe_select(fens,fes,struct ('label', 13))
%
% flood
%
% Select all FEs connected together (Starting from node 13):
%     fe_select(fens,fes,struct ('flood', true, 'startfen', 13))
%
% facing
%
% Select all FEs "facing" in the direction [1,0,0] (along the x-axis):
%     fe_select(fens,fes,struct ('facing', true, 'direction', [1,0,0]))
% Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
%     fe_select(fens,fes,struct ('facing', true, 'direction', @(x)(+[x(1:2),0])))
% Here x is the centroid of the nodes of each selected fe.
% Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
%     fe_select(fens,fes,struct ('facing',1,'direction',@(x)(+[x(1:2),0]),'tolerance',0.01))
% Here the fe is considered facing in the given direction if the dot
% product of its normal and the direction vector is greater than tolerance.
%
% The option 'inflate' may be used to increase or decrease the extent of
% the box (or the distance) to make sure some nodes which would be on the
% boundary are either excluded or included.
%
% nearestto
%
% Select the geometric cell nearest to the given location.
% nearest=fe_select(fens,fes,struct('nearestto',[2.3,-2]))
%
% distance
%
% Select the finite elements within a given distance from the center.
% Example: fe_select(fens,fes,struct ('distance',0.5, 'from',[1 -1])), selects
% elements whose nodes are Less than 0.5 units removed from the point [1 -1].
%
% smoothpatch
%
% Select all FEs that are part of a smooth surface. 
% For instance, starting from the finite element number 13, select all finite elements 
% whose normals defer from the normal of the Neighbor element by less than
% 0.05 in the sense that dot(n1,n2)>sqrt(1-0.05^2).
%     fe_select(fens,fes,struct ('flood', true, 'startfe', 13,  'normaldelta',0.05))
% % Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
%     fe_select(fens,fes,struct ('facing', true, 'direction', @(x)(+[x(1:2),0])))
% Here x is the centroid of the nodes of each selected fe.
% % Select all FEs "facing" in the direction [x(1),x(2),0] (away from the z-axis):
%     fe_select(fens,fes,struct ('facing',1,'direction',@(x)(+[x(1:2),0]),'tolerance',0.01))
% Here the fe is considered facing in the given direction if the dot
% product of its normal and the direction vector is greater than tolerance.
%
%
% Output:
% felist= list of finite elements from the set that satisfy the criteria
%
% Examples:
%     topl =fe_select (fens,bfes,struct('box', [0,L,-Inf,Inf,b/2,b/2],...
%                             'inflate',tolerance));
%     %  The subset  of the finite elements that are in the box is
%     traction.fes= subset(bfes,topl);
%
% Note: This function uses the node selection function fenode_select() to search
% the nodes.
%
% See also: fenode_select
%
inside = true;
anynode = false;
flood = false;
facing = false;
smoothpatch= false;
label = [];
overlapping_box = false;
nearestto = false;
if isfield(options,'anynode')
    inside = false;
    anynode = true;
end
if isfield(options,'flood')
    flood = true;
    if isfield(options,'startfen')
        startfen = options.startfen;
    else
        error('Need the identifier of the Starting node, startfen');
    end
end
if isfield(options,'facing')
    facing = true;
    if isfield(options,'direction')
        direction = options.direction;
    else
        error('Need the direction as a three element array');
    end
    tolerance =0;
    if isfield(options,'tolerance')
        tolerance = options.tolerance;
    end
end
if isfield(options,'smoothpatch')
    smoothpatch = true;
    if isfield(options,'normaldelta')
        normaldelta = options.normaldelta;
    else
        normaldelta= 0.05;
    end
    if isfield(options,'startfe')
        startfe = options.startfe;
    else
        error('Need the number of the starting finite element');
    end
end
if isfield(options,'label')
    label = options.label;
end
if isfield(options,'overlapping_box')
    overlapping_box = true;
    box=options.overlapping_box;
    bounding_boxes =[];% precomputed bounding boxes of fes can be passed in
    if isfield(options,'bounding_boxes')
        bounding_boxes = options.bounding_boxes;
    end
    inflate =0;
    if isfield(options,'inflate')
        inflate = (options.inflate);
    end
    box = inflate_box (box, inflate);
end
if isfield(options, 'nearestto')
    nearestto = true;
    locations = options.nearestto;
end

%     Select based on fe label
if ~isempty(label)
    lbl=fes.label;
    felist= [];
    if (isempty(lbl))
        return;
    end
    conn=fes.conn;
    for i=1:size(conn,1)
        if label==lbl(i)
            felist(end+1) =i;
        end
    end
    return;
end

% Select by flooding
if (flood)
    fen2fe_map=fenode_to_fe_map (struct ('fes',fes));
    gmap=fen2fe_map.map;
    conn=fes.conn;
    felist= zeros(1, size(conn,1));
    felist(gmap{startfen})=1;
    while true
        pfelist=felist;
        markedl=find(felist~=0);
        for j=markedl
            for k=conn(j,:)
                felist(gmap{k})=1;
            end
        end
        if sum(pfelist-felist)==0, break; end
    end
    felist =find(felist~=0);
    return;
end

% Select by in which direction the normal of the fes face
if (facing)
    felist= [];
    xs =fens.xyz;
    sdim =size(xs,2);
    mdim=fes.dim;
    if (mdim~=sdim-1)
        error ('"Facing" selection of fes make sense only for Manifold dimension ==Space dimension-1')
    end
    param_coords =zeros(1,mdim);
    Need_Evaluation = (strcmp(class (direction),'inline') || strcmp(class (direction),'function_handle'));
    if ( ~Need_Evaluation)
        d = reshape(direction,1,[])/norm(direction);
    end
    Nder = bfundpar (fes, param_coords);
    conns=fes.conn;
    for i=1: size(conns,1)
        conn=conns(i,:);
        xyz =zeros(length (conn),sdim);
        xyz =xs(conn,:);
        Tangents =xyz'*Nder;
        N = normal (Tangents,sdim, mdim);
        if (Need_Evaluation)
            d=feval(direction,mean(xyz));
            d = reshape(d,1,[])/norm(d);
        end
        if (dot(N,d)>tolerance)
            felist(end+1)=i;
        end
    end
    felist =unique(felist);
    return;
end


% Select by the change in normal
if (smoothpatch)
    mincos=sqrt(1-normaldelta^2);
    xs=fens.xyz;
    sdim =size(xs,2);
    mdim=fes.dim;
    if (mdim~=sdim-1)
        error ('"Smoothpatch" selection of fes make sense only for Manifold dimension ==Space dimension-1')
    end
    fen2fe_map=fenode_to_fe_map (struct ('fes',fes));
    gmap=fen2fe_map.map;
    conn=fes.conn;
    param_coords =zeros(1,mdim);% This is a hack: this may not be the proper location for all element types
    Nder = bfundpar (fes, param_coords);
    felist= zeros(1, size(conn,1));
    felist(startfe)=1;
    while true
        pfelist=felist;
        markedl=find(felist~=0);
        for j=markedl
            xyz =xs(conn(j,:),:);
            Tangents =xyz'*Nder;
            Nj = normal (Tangents,sdim, mdim);
            for k=conn(j,:)
                for ke=gmap{k}
                    xyz =xs(conn(ke,:),:);
                    Tangents =xyz'*Nder;
                    Nk = normal (Tangents,sdim, mdim);
                    if (dot(Nj,Nk)>mincos)
                        felist(ke)=1;
                    end
                end
            end
        end
        if sum(pfelist-felist)==0, break; end
    end
    felist =find(felist~=0);
    return;
end
    

% Select all FEs whose bounding box overlaps given box
if (overlapping_box)
    felist= [];
    xs =fens.xyz;
    conns=fes.conn;
    if (isempty(bounding_boxes))
        bounding_boxes =zeros(size(conns,1),length(box));
        for i=1: size(conns,1)
            conn=conns(i,:);
            xyz =xs(conn,:);
            bounding_boxes(i,:)= bounding_box(xyz);
        end
    end
    for i=1: size(conns,1)
        if (boxes_overlap (box,bounding_boxes(i,:)))
            felist(end+1)=i;
        end
    end
    felist =unique(felist);
    return;
end


% get the fe nearest to the supplied point
if (nearestto)
    
    felist = [];
    
    % get the locations of all the nodes
    xs = fens.xyz;
    %     Disqualify all the nodes that are not connected to the finite
    %     elements on input 
    cn = connected_nodes(fes);
    xs(setdiff(1:size(xs,1),cn),:) =Inf;
    
    % get the connectivity of all the FEs
    conns = fes.conn;
    
    for i = 1:size(locations,1)
        
        % Get the smallest distances between the node locations and the
        % desired location using ipdm (by John D'Errico)
        distance = ipdm(xs, locations(i,:));
        
        % Get array index of smallest distance to location from
        % all nodes
        [junk,IX] = min(distance);
        
        % Find the fes connected to the nearest node
        [crows, ccols] = find(conns == IX(1));
        
        if numel(crows) == 1
            % if there's only one, this is superb, we can add it to the
            % list and continue to the next location
            felist(end+1) = crows;
        else
            
            % otherwise we must determine which cell is closest based on
            % the other connected nodes
            nearconns = conns(crows, :);
            
            nearconndist = zeros(size(nearconns));
            
            % Get the distances between all the nodes in the nearest
            % FEs and the nearest node to the location
            for j = 1:size(nearconns, 1)
                for k = 1:size(nearconns, 2)
                    nearconndist(j,k) = ipdm(xs(nearconns(j,k), :), locations(i,:));
                end
            end
            
            % set the distance from the nearest node to in the cell to
            % the location to infinity
            nearconndist(nearconns == IX(1)) = inf;
            
            % order the FEs by the smallest distance of any connected
            % node from the nearest node
            [junk,IX] = sort(min(nearconndist,[],2));
            
            % return the index in conns to the first of these ordered
            % FEs
            if (~isempty(crows))
                felist(end+1) = crows(IX(1));
            end
            
        end
        
    end
    
    return;
    
end

%     Select based on location of nodes
nodelist=fenode_select(fens, options);
conn=fes.conn;
felist= [];
for i=1: size(conn,1)
    tf = ismember(conn(i,:), nodelist);
    if inside
        if sum(tf) == length(conn(i,:))
            felist(end+1) =i;
        end
    end
    if anynode
        if sum(tf) >= 1
            felist(end+1) =i;
        end
    end
end

end

function N = normal (Tangents,sdim, mdim)
[sdim, mdim] = size(Tangents);
if     mdim==1 % 1-D fe
    N = [Tangents(2,1),-Tangents(1,1)];
elseif     mdim==2 % 2-D fe
    N =skewmat(Tangents(:,1))*Tangents(:,2);
else
    error('Got an incorrect size of tangents');
end
N=N/norm(N);
end

