function map = t4util_v_to_t4e_map (nv,e)
% Construct map from vertices to edges.
%
% function map = t4util_v_to_t4e_map (nv,e)
%
    map = cell(nv,1);
    for i=1:length( map)
        map{i}=[];
    end
    for i=1:size(e,1)
        conn =e(i,:);
        for j=1:length(conn)
            m=map{conn(j)};
            m=[m i];
            map{conn(j)}=m;
        end
    end
end