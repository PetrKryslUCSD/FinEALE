% Make spherical surfaces
function test_spherical_surface1
axis equal
Radius=300;
nperradius=4;
[xyz,conn]=make_spherical_surface_t3(Radius,nperradius);
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');

Radius=300*2;
nperradius=4;
[xyz,conn]=make_spherical_surface_t3(Radius,nperradius);
xyz= xyz +ones(size( xyz,1 ),1)*[1.2*Radius,-Radius,Radius];
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'b','EdgeColor','k');
end