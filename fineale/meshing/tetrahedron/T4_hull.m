function [fens,fes] = T4_hull(vertices)
    % Tetrahedral (T4) mesh of a cloud of points.
    %
    % function [fens,fes] = T4_hull(vertices)
    %
    % Arguments
    % vertices = one row per vertex, coordinates of vertices (points).
    %
    % The function generates a tetrahedral mesh that is the convex hull
    % of the vertices supplied on input.
    %
    %
    % Examples: 
    %     [fens,fes] = T4_cylinderdel(3.1,1.7,5,2);
    %     figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
    %     [fens,fes] =  T4_hull(fens.xyz);
    %     figure; drawmesh({fens,fes},'fes','facecolor','b'); hold on

    
    T = delaunayn(vertices);
    fes=fe_set_T4(struct ('conn',T));
    fens=fenode_set (struct('xyz',vertices));
end

%compute volume of a tetrahedron
% Given the 4x3 vertex coordinate matrix V of a tetrahedron, TETVOL(V)
% returns the volume of the tetrahedron.
function vol = tetvol(v)
    vol = det([v(2,:)-v(1,:);v(3,:)-v(1,:);v(4,:)-v(1,:)])/6;
    return;
end
