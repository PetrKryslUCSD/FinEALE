% Compute tessellation of an isosurface within hexahedral finite elements.
%
% function  [f,v]= isosurf_faces(feSet, conns, geom, u, scalarfield, isovalue, nref)
% 
%    feSet= finite element set object for which the algorithm is getting called
%    conns= array of connectivity of hexahedral  shapes, eight columns 
%           and as many rows as needed 
%    geom= geometry field
%    u = displacement field 
%    scalarfield =field with scalar vertex data,
%    isovalue = value of the isosurface
%    color = what color should the isosurface be?  Default: [1,1,1]/2;
%    nref = how refined should the parametric space of the geometric cell
%         be? Default: 2.
%
function  [f,v]= isosurf_faces(feSet, conns, geom, u, scalarfield, isovalue, nref)
    f= []; v= []; numv=0;
    maxconn=max(max(conns));% largest node number
    ss =gather_values(scalarfield, (1:maxconn));% values of the scalar field
    if ((min(ss)>isovalue) ||  (max(ss)<isovalue))
        return% nothing to be done
    end
    
    xis =linspace(-1,+1,nref);
    [X,Y,Z] = ndgrid(xis,xis,xis);
    
    Nstransp={};
    for k=1:length(xis)
        for j=1:length(xis)
            for i=1:length(xis)
                Nstransp{i, j, k} = (bfun(feSet,  [xis(i),xis(j),xis(k)]))';
            end
        end
    end
    
    xs = gather_values(geom, (1:maxconn)); % coordinates of nodes
    us = gather_values(u, (1:maxconn)); % displacements of nodes
    xus=xs+us;
    
    for jconn=1:size(conns,1)
        conn = conns(jconn,:); % connectivity
        scalars=ss(conn);
        if ~((min(scalars)>isovalue) ||  (max(scalars)<isovalue))
            V=0*Z;
            for k=1:length(xis)
                for j=1:length(xis)
                    for i=1:length(xis)
                        V(i,j,k)=Nstransp{i, j, k}*scalars;
                    end
                end
            end
            [nf,nv] = isosurface(X,Y,Z,V,isovalue);
            xu=xus(conn,:);
            for i=1:size(nv,1)
                N = bfun(feSet,  nv(i,:));
                nv(i,:)=N'*xu;
            end
            numv=size(v,1);
            f= [f;nf+numv];
            v= [v;nv];
        end % there is an iso-surface within this volume
    end % loop over all volumes
end % done