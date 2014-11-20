% Make one half of a hemisphere.
function mkhemisphere
    R= 300;
    Ri= 1;
    thickness =R-Ri;
    nlayers =13;
    [fens,fes]=Q4_sphere(R,1,thickness);
    function xyz= radially(xyz, layer)
        xyz= (Ri+layer/nlayers*thickness)*xyz/norm(xyz);
    end
    [fens,fes] = H8_extrude_Q4(fens,fes,nlayers,@radially);
    
    [fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0]);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, Ri/1000);
    fes=cat(fes1,fes2);
    drawmesh({fens,fes},'fes','facecolor','red')
    %         function xyz= up(xyz, layer)
    %         xyz= [xyz+(layer^2)*[0.2,-0.2], layer*2.5];
    %     end
    %     [fens,fes] = extrudeq4(fens,fes,3,@up);
    %     % drawmesh({fens,fes},'fes','nodes','facecolor','red')
    %     drawmesh({fens,fes},'fes','facecolor','red')
    % ix =gcell_select(fens,bg,...
    %     struct ('box',[-100 100 -100 0 -100 0],'inflate', 0.5))
    % drawmesh({fens,bg(ix)},'facecolor','red','shrink', 1.0)
end