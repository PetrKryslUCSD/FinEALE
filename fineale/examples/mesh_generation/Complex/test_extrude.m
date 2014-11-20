%Test extrusion of quadrilaterals into hexahedra.
function test_extrude
    [fens,fes] = Q4_L2x2;
    function xyz= up(xyz, layer)
        xyz= [xyz+(layer^2)*[0.2,-0.2], layer*2.5];
    end
    [fens,fes] = H8_extrude_Q4(fens,fes,3,@up);
    % drawmesh({fens,fes},'fes','nodes','facecolor','red')
    drawmesh({fens,fes},'fes','facecolor','red')
    % ix =fe_select(fens,bg,...
    %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
    % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
end