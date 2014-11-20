function map = t4util_v_to_t4_map (nv,t)
% Construct map from vertices to tetrahedra.
%
% function map = t4util_v_to_t4_map (nv,t)
%
    map = cell(nv,1);
    for i=1:length( map)
        map{i}=[];
    end
    for i=1:size(t,1)
        conn =t(i,:);
        for j=1:length(conn)
            m=map{conn(j)};
            m=[m i];
            map{conn(j)}=m;
        end
    end
end