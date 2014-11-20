% Test tetrahedral blocks.
function test_plate_surface1
axis equal
Length=500;
Width=Length/2;
Height=Length/4;
nperLength=5;
nperWidth=3;
nperHeight=2;
[fens,fes] = T4_blocka(Length,Width,Height,nperLength,nperWidth,nperHeight);
fes=mesh_boundary(fes);
xyz=fens.xyz;
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');


[fens,fes] = T4_blocka(Length,Width,Height,nperLength,nperWidth,nperHeight);
fes=mesh_boundary(fes);
xyz=fens.xyz;
xyz= xyz*rotmat([pi/4,0,0])';
xyz= xyz+ones(size( xyz,1 ),1)*[Length,-Length,1.2*Length];
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'b','EdgeColor','k');
view (3)
end