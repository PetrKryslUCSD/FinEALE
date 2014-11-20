function f = t4util_interfaces(t,tmid,tmid12)
% Extract the interface triangles from a T4 connectivity.
%
% f = t4util_interfaces(t,tmid,tmid12)
%
% t= connectivity of the tetrahedra
% tmid= material identifier of each tetrahedron
% tmid12= array of two material identifiers between which the interface
% exists
%
% The function extracts the connectivity of the interface triangles in
% multi-material meshes. 
% 
    f= [];
    for j=1:size(tmid12,1)
        tmid1 =tmid12(j,1);
        tmid2 =tmid12(j,2);
        f1 = t4util_bdry(t(find(tmid==tmid1),:));
        f2 = t4util_bdry(t(find(tmid==tmid2),:));
        [c, ia, ib]  = intersect(sort(f1,2),sort(f2,2),'rows');
        f =[f;f1(ia,:)];
    end
end
