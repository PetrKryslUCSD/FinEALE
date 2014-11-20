% Test extrusion of lines into quadrilaterals
function test_extrudel2
    [fens,fes] = L2_block(4.0,5, struct('other_dimension', 3.1));
    function xyz= upz(xyz, layer)
        xyz= [xyz+(layer^2)*[0.,-0.2], layer*2.5];
    end
    function xyz= up(xyz, layer)
        xyz= [xyz(1),layer*.5];
    end
    [fens,fes] = Q4_extrude_L2(fens,fes,2,@up)
    % drawmesh({fens,fes},'fes','nodes','facecolor','red')
    drawmesh({fens,fes},'fes','facecolor','red')
    % ix =fe_select(fens,bg,...
    %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
    % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
end