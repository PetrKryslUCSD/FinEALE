% Test merging of meshes.
function test_merging
    
    L=2;H=3;
    
    Result=true;
    
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{2},fesa{2}] = Q4_block(L,H,3, 4, 1.0);
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.);
    Result=Result&&(count(fens)==40);
    
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{2},fesa{2}] = Q4_block(L,H,3, 4, 1.0);
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.01);
    Result=Result&&(count(fens)==20);
      
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{2},fesa{2}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{3},fesa{3}] = Q4_block(L,H,3, 4, 1.0);
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.);
    Result=Result&&(count(fens)==60);
    
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{2},fesa{2}] = Q4_block(L,H,3, 4, 1.0);
    [fensa{3},fesa{3}] = Q4_block(L,H,3, 4, 1.0);
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.01);
    Result=Result&&(count(fens)==20);
    
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    fensa{2} = transform_apply(fensa{1},@(x,d)x+[L,0], []);
    fesa{2} =fesa{1};
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.01);
    Result=Result&&(count(fens)==40-5);
    
    
    [fensa{1},fesa{1}] = Q4_block(L,H,3, 4, 1.0);
    fensa{2} = transform_apply(fensa{1},@(x,d)x+[L,0], []);
    fesa{2} =fesa{1};
    fensa{3} = transform_apply(fensa{1},@(x,d)x+[0,H], []);
    fesa{3} =fesa{1};
    [fens,fesa] = merge_n_meshes(fensa, fesa, 0.01);
    Result=Result&&(count(fens)==60-5-4);
    
    gv=drawmesh({fens,fesa{1}},'nodes'); 
    for j=1:length(fesa)
    gv=drawmesh({fens,fesa{j}},'gv',gv,'fes'); 
    end 
    
    Result
    assignin('caller','fineale_test_passed',Result);
end