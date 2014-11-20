% Create triangular meshes on elliptical cross-section cylinders.
function test_elliptcylinder_surface1
axis equal
Length=500;
MinorRadius=Length/4;
Major_to_minor = 2.3;
nperLength=5;
nperRadius=4;
[fens,fes] = H8_cylinder_n(MinorRadius, Length, nperRadius, nperLength);
fes=mesh_boundary(fes);
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
xyz=fens.xyz;
xyz(:,1)=Major_to_minor*xyz(:,1);
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');

[fens,fes] = H8_cylinder_n(MinorRadius, Length, nperRadius, nperLength);
fes=mesh_boundary(fes);
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
xyz=fens.xyz;
xyz(:,1)=Major_to_minor*xyz(:,1);
xyz= xyz*rotmat([0,0,pi/2])';
xyz= xyz+ones(size( xyz,1 ),1)*[Length,-Length,1.2*Length];
conn=fes.conn;
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'b','EdgeColor','k');
view(3)
end