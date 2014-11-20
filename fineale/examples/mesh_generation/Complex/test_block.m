% Calculate the visual table comparing the different block-generation methods..
function test_block
    L=2;H=3;
    Handles ={@Q4_block,@T3_block,@T3_blockc,@T3_ablock,@T3_ablockc,@T3_cblock,@T3_crossblock,@T3_blockn,@T3_ablockn,@T3_blocku};
    nr=2;
    nc =ceil(length( Handles )/2);
    gv= graphic_viewer;
    i=1;
    for m= Handles
        subplot(nr,nc,i); i=i+1;
         axis equal ;
        [fens,fes] = m{1}(L,H,3, 4, 1.0);
        drawmesh({fens,fes},'gv',gv,'fes','facecolor','red');
        title(strrep(func2str (m{1}),'_', '\_'))
        view(2)
        pause(2)
    end
end