% Test generation of a general hexahedron.
xyz = [3, 1, 6; -5, 2, 1];
[fens,fes] = H8_hexahedron(xyz,12,3,4);
%     [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens2, fes2, eps);
%     fes= cat(2,fes1,fes2);
    % drawmesh({fens,fes},'fes','nodes','facecolor','red')
    drawmesh({fens,fes},'fes','facecolor','red'); hold on
%     drawmesh({fens,fes2},'fes','facecolor','Green'); hold on
    % ix =fe_select(fens,bg,...
    %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
    % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
