% Compute the line segments on level curves of a nodal field.
%
% function  [e,v]= isosurf_edges(conns, geom, u, scalarfield, isovalue)
% 
%    conns= array of connectivity of triangles, three columns and as many rows as needed
%    geom= geometry field
%    u = displacement field 
%    scalarfield =field with scalar vertex data,
%    isovalue = value of the isosurface
%    color = what color should the isosurface be?  Default: [1,1,1]/2;
%
function  [e,v]= isosurf_edges(conns, geom, u, scalarfield, isovalue)
e= []; v= []; numv=0;
maxconn=max(max(conns));% largest node number
ss =gather_values(scalarfield, (1:maxconn));% values of the scalar field
if ((min(ss)>isovalue) ||  (max(ss)<isovalue))
    return% nothing to be done
end



xs = gather_values(geom, (1:maxconn)); % coordinates of nodes
us = gather_values(u, (1:maxconn)); % displacements of nodes
xus=xs+us;
nv = zeros(2,size(xus,2));
for jconn=1:size(conns,1)
    conn = conns(jconn,:); % connectivity
    scalars=ss(conn);
    if ~((min(scalars)>isovalue) ||  (max(scalars)<isovalue))
        
        Pcont = tricontour_edge_([0,0;1,0;0,1],[1,2,3],scalars,isovalue);
        xu=xus(conn,:);
        
        N = [1-Pcont(1,1)-Pcont(1,2),Pcont(1,1),Pcont(1,2)];
        nv(1,:)=N*xu;
        N = [1-Pcont(2,1)-Pcont(2,2),Pcont(2,1),Pcont(2,2)];
        nv(2,:)=N*xu;
        
        numv=size(v,1);
        e= [e;[1,2]+numv];
        v= [v;nv];
    end % there is an iso-surface within this volume
end % loop over all volumes
end % done


function Pcont = tricontour_edge_(p,t,Hn,T)
Pcont = []; %all points [x y] of the contour

    Pt = []; %all points in the triangle
    e = [t([1,2]); t([2,3]); t([3,1])]; %edges in each triangle
    
    %looping through the edges in the m'th triangle
    for k = 1:3
        
        %determine whether or not the contour will pass through the
        %edge, will go through all 3 edges
        %If the contour passes through the node, it will ==, probably
        %want to change something.
        if T >= min(Hn(e(k,:))) && T <= max(Hn(e(k,:)))
            
            r = (T - Hn(e(k,1)))/(Hn(e(k,2))-Hn(e(k,1)));
            %r might be NaN if both nodes are ==, 0/0
            
            P = p(e(k,1),:) + r*(p(e(k,2),:) - p(e(k,1),:));
            
            Pt = [Pt;P];
            
        end
        
        %check for special case where contour passes through a node,
        %will delete repeating point.
        %Probably can combine this above while checking if there is a
        %point on each edge when node==T(i)
        if size(Pt,1) == 3
            if min(Pt(1,:) == Pt(2,:));
                Pt(1,:) = [];
            elseif min(Pt(2,:) == Pt(3,:));
                Pt(2,:) = [];
            else min(Pt(1,:) == Pt(3,:));
                Pt(1,:) = [];
            end
        end
        
    end
    
    Pcont = [Pcont;Pt];
    
end