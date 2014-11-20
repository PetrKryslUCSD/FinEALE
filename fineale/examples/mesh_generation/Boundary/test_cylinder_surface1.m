% Create meshes of cylindrical surfaces.
function test_cylinder_surface1
axis equal
Length=500;
Radius=Length/4;
nperLength=5;
nperRadius=4;
[fens,fes] = H8_cylinder_n(Radius, Length, nperRadius, nperLength);
fes=mesh_boundary(fes);
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
xyz=fens.xyz;
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');

[fens,fes] = H8_cylinder_n(Radius, Length, nperRadius, nperLength);
fes=mesh_boundary(fes);
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
xyz=fens.xyz;
xyz= xyz*rotmat([pi/4,0,0])';
xyz= xyz+ones(size( xyz,1 ),1)*[Length,-Length,1.2*Length];
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'b','EdgeColor','k');

[fens,fes] = H8_cylinder_n(Radius, Length, nperRadius, nperLength);
fes=mesh_boundary(fes);
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
xyz=fens.xyz;
xyz= xyz+ones(size( xyz,1 ),1)*[0,0,-Length/2];
xyz= xyz*rotmat([0,pi/4,0])';
xyz= xyz+ones(size( xyz,1 ),1)*[Length,0,0];
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'g','EdgeColor','k');
end