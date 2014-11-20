function vn =  vertex_neighbors(vn,f,v)
% Find the node neighbors in the mesh.
%
% function vn =  vertex_neighbors(vn,f,v)
%
% vn= cell array, element I holds an array of numbers of nodes  
%     which are connected to node I (including node I).  When this array is 
%     supplied as input the information from the current call is added to 
%     the array vn; otherwise (when vn is empty on input) the array is created 
%     and returned.
% f= connectivity of the mesh, one row per element
% v= locations of the nodes, three columns, one row per node
    if (isempty(vn))
        vn=cell(1,size(v,1));
    end
    for I= 1:size(f,1)
        for r= 1:size(f,2)
            vn{f(I,r)}=[vn{f(I,r)},f(I,:)];
        end
    end
    for I= 1:size(vn,2)
        vn{I}=unique(vn{I});
    end
end