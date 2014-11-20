function [t,v] =smoothing(t,v,options)
% General smoothing of meshes.
%
% function [t,v] =smoothing(t,v,options)
%
% Fields of the structure options, all are optional:
% method='laplace' or 'taubin' (Default is 'taubin'.)
% f=boundary faces (optional)
% bv=boundary vertices (optional)
% bv_from_f=compute boundary of vertices from the boundary faces, true or
% false.  Tetrahedra and hexahedra are supported.
% npass=how many passes of smoothing? default is 2.
    iv=v;
    fvn = [];
    bv = [];
    bv_from_f = false;
    f= [];
    npass =2;
    method =@taubin_smoother;
    if (exist('options','var'))
        if (isfield(options,'bv'))
            bv= options.bv;
        end
        if (isfield(options,'f'))
            f= options.f;
        end
        if (isfield(options,'bv_from_f'))
            bv_from_f= options.bv_from_f;
        end
        if (isfield(options,'npass'))
            npass= options.npass;
        end
        if (isfield(options,'method'))
            if ( strcmp(options.  method,'laplace'))
                method =@laplace_smoother;
            elseif ( strcmp(options.  method,'taubin'))
                method =@taubin_smoother;
            end
        end

    end

    if (isempty(bv))
        bv=0*v(:,1);
    end
    if (bv_from_f)
        if (isempty(f))
            if ( size(t,2)==4)
                f = t4util_bdry(t); % Extract the boundary
            elseif ( size(t,2)==8)
                f = h8util_bdry(t); % Extract the boundary
            else
                f=[];
            end
            bv(f)=1;
        end
    end
   
    % find neighbors for the Surface connections
    fvn = vertex_neighbors([],f,v);
    % Smoothing considering only surface connections
    fv =  method(v,fvn,0*bv,npass,0.5,-0.5);
    
    % find neighbors for the Volume connections
    vn =  vertex_neighbors([],t,v);
    % Smoothing considering all connections through the volume
    v =  method(v,vn,bv,npass,0.5,-0.5);

    % Correction of the vertices of the surface
    if (~isempty(fvn))
        for i= 1:length(fvn)
            if (~isempty(fvn {i}))
                v(i,:)= fv(i,:);
            end
        end
    end
end

function nv =  taubin_smoother(v,vn,bv,npass,lambda,mu)
    nv=v;
    for I= 1:npass
        o=randperm(length(vn));
        damping_factor=lambda;
        for k= 1:length(vn)
            r=o(k);
            n=vn{r};
            if (~isempty(n)) && (length(n)>1) && (~bv(r))
                nv(r,:)=(1-damping_factor)*v(r,:)+ damping_factor*(sum(v(n,:))-v(r,:))/(length(n)-1);
            end
        end
        v=nv;
        damping_factor=mu;
        for k= 1:length(vn)
            r=o(k);
            n=vn{r};
            if (~isempty(n)) && (length(n)>1) && (~bv(r))
                nv(r,:)=(1-damping_factor)*v(r,:)+ damping_factor*(sum(v(n,:))-v(r,:))/(length(n)-1);
            end
        end
        v=nv;
    end
end

function nv =  laplace_smoother(v,vn,bv,npass,lambda,mu)
    damping_factor=lambda;
    for I= 1:npass
        nv=v;
        o=randperm(length(vn));
        for k= 1:length(vn)
            r=o(k);
            n=vn{r};
            if (~isempty(n)) && (length(n)>1) && (~bv(r))
                nv(r,:)=(1-damping_factor)*v(r,:)+ damping_factor*(sum(v(n,:))-v(r,:))/(length(n)-1);
            end
        end
        v=nv;
    end

end



