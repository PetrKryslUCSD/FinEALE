% Construct arrays to describe a tetrahedron mesh created from voxel image.
%
% function [t,v,tmid] = voximg_to_T4(img,voxval)
%
% img = 3-D image (array), the voxel values  are arbitrary
% voxval =range of voxel values to be included in the mesh, voxval =
% [minimum value, maximum value].  Minimum value == maximum value is
% allowed.
% Output:
% t = array of tetrahedral connectivities, one tetrahedron per row
% v =Array of vertex locations, one vertex per row
function [t,v,tmid] = voximg_to_T4(img,voxval)
    M=size(img, 1); N=size(img, 2); P=size(img, 3);
    Nvoxval =length(find((img>=min(voxval))  & (img<=max(voxval))));
    t =zeros(5*Nvoxval,4,'int32');
    v =zeros(prod((size(img)+1)),3,'int16');
    tmid =zeros(5*Nvoxval,1,'int16');
    t4ia = [8, 4, 7, 5; 6, 7, 2, 5; 3, 4, 2, 7; 1, 2, 4, 5; 7, 4, 2, 5];
    t4ib = [7, 3, 6, 8; 5, 8, 6, 1; 2, 3, 1, 6; 4, 1, 3, 8; 6, 3, 1, 8];
    
    Slice =zeros(2,N+1,P+1);
    nv =0;
    nt =0;
    for I= 1:M
        tf =ismember(img(I,:,:), voxval);
        for J= 1:N
            for K= 1:P
                if (tf(1,J,K))
                    store_hex (I,J,K);
                end
            end
        end
        Slice(1,:,:) =Slice(2,:,:) ;
        Slice(2,:,:) =0;
    end
    v=v(1:nv,:);
    t=t(1:nt,:) ;
    tmid=tmid(1:nt) ;
    v=double(v)-1;% Adjust the locations so that the corner of the volume is at 0,0,0
    return;

    function vidx = find_vertex (I,IJK)
        vidx = zeros(1,size(IJK,1));
        for r= 1:size(IJK,1)
            if (Slice(IJK(r,1),IJK(r,2),IJK(r,3))==0)
                nv=nv+1;
                v(nv,:) =IJK(r,:)+ [I-1,0,0];
                Slice(IJK(r,1),IJK(r,2),IJK(r,3)) =nv;
            end
            vidx(r) =Slice(IJK(r,1),IJK(r,2),IJK(r,3));
        end
    end


    function store_hex (I,J,K)
%         locs =[I,J,K;I+1,J,K;I+1,J+1,K;I,J+1,K];[I,J,K+1;I+1,J,K+1;I+1,J+1,K+1;I,J+1,K+1];
        locs =[1,J,K;1+1,J,K;1+1,J+1,K;1,J+1,K;1,J,K+1;1+1,J,K+1;1+1,J+1,K+1;1,J+1,K+1];
        vidx = find_vertex (I,locs);
        for r=1:5
            nt =nt +1;
            if (mod (sum( [I,J,K] ),2)==0)
                t(nt,:) =vidx(t4ia(r,:));
            else
                t(nt,:) =vidx(t4ib(r,:));
            end
            tmid(nt) =img(I,J,K);
        end
    end
end