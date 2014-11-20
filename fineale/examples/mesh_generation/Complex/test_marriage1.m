% Test merging of extruded meshes.
function test_marriage1
    [fens,fes] = Q4_L2x2;
    function xyz= up(xyz, layer)
        xyz= [xyz+(layer^2)*[0.2,-0.2], layer*2.5];
    end
    function xyz= down(xyz, layer)
        xyz= [xyz-(layer^2)*[0.2,-0.2], -layer*2.5];
    end
    [fens1,fes1] = H8_extrude_Q4(fens,fes,4,@up);
    [fens2,fes2] = H8_extrude_Q4(fens,fes,5,@down);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens2, fes2, eps);
    fes= cat(fes1,fes2);
    % drawmesh({fens,fes},'fes','nodes','facecolor','red')
    drawmesh({fens,fes},'fes','facecolor','red'); hold on
%     drawmesh({fens,fes2},'fes','facecolor','Green'); hold on
    % ix =fe_select(fens,bg,...
    %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
    % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
    
end