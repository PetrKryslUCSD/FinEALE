% Construct meshes for the elliptical hole configuration.
function test_elliphole
    xradius=1;yradius=2;L=5;H=8;nL=5;nH=3;nR=5;
    options.thickness =1.0;
    [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options);
    drawmesh({fens,fes},'fes','facecolor','red');
    view(2)
    pause(2)
xradius=2;yradius=2;L=5;H=8;nL=5;nH=7;nR=9;
    options.thickness =1.0;
    [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options);
    drawmesh({fens,fes},'fes','facecolor','red');
    view(2)
    pause(2)
xradius=2;yradius=2;L=5;H=8;nL=1;nH=1;nR=1;
    options.thickness =1.0;
    [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options);
    drawmesh({fens,fes},'fes','facecolor','red');
    view(2)
    pause(2)
