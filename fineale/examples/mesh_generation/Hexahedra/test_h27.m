% Test generation of a 27-node hexahedral block.
[fens,fes] = H27_block(12,3,4,1,1,1);
fens.xyz= [-1,-1,-1;
                    1,-1,-1;
                    1,1,-1;
                    -1,1,-1;-1,-1,1;
                    1,-1,1;
                    1,1,1;
                    -1,1,1;0,-1,-1;
                    1,0,-1;
                    0,1,-1;
                    -1,0,-1;0,-1,1;
                    1,0,1;
                    0,1,1;
                    -1,0,1;-1,-1,0;
                    1,-1,0;
                    1,1,0;
                    -1,1,0;0,0,-1;
                    0,-1,0;
                    1,0,0;
                    0,1,0;
                    -1,0,0;
                    0,0,1;0,0,0];
                fes.conn=1:27;
gv=drawmesh({fens,fes},'nodes','fes','facecolor','none','offset',0.04); hold on
bg=mesh_boundary(fes);
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [-1,0,0],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y'); hold on
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [0,-1,0],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y'); hold on
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,-1,],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','y'); hold on
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [+1,0,0],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','r'); hold on
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [0,+1,0],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','g'); hold on
l=[fe_select(fens,bg,struct ('facing', true, 'direction', [0,0,+1,],'tolerance',0.01))];
gv=drawmesh({fens,subset(bg,l)},'gv',gv,'fes','facecolor','b'); hold on
labels