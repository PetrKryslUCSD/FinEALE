% Create surface mesh of a cylinder capped with hemispheres.
function test_cappedcylinder_surface1
axis equal
Length=500;
Radius=Length/4;
nperLength=8;
nperRadius=6;
tol=Radius/nperRadius/2/100;
% Make 1/8 of the spherical surface
[fens,fes]=Q4_sphere_n(Radius,nperRadius,1);
%     Mirror to get 1/4
[fens1,fes1] = mirror_mesh(fens, fes, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
fes=cat(fes1,fes2);
%     Mirror again to get 1/2: this is the top cap
[fens1,fes1] = mirror_mesh(fens, fes, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tol);
tfes=cat(fes1,fes2);
tfens=fens;
%     Mirror to get the bottom cap
[bfens,bfes] = mirror_mesh(tfens, tfes, [0,0,-1], [0,0,0], @(c)c([1, 4, 3, 2]));
% Extract the circle boundary
bcfes=mesh_boundary(bfes);
% extrude the circular boundary to create the cylindrical surface
[cfens,cfes] = Q4_extrude_L2(fens,bcfes,nperLength,...
    @(xyz, layer)xyz+[0,0,(layer/nperLength)*Length]);
% Translate the top cap
tfens = transform_apply(tfens,@(x,d)x+[0,0,Length], []);
% Merge everything together
[fens,fes1,fes2] = merge_meshes(tfens, tfes, bfens, bfes, tol);
fes=cat(fes1,fes2);
[fens,fes1,fes2] = merge_meshes(fens, fes, cfens, cfes, tol);
fes=cat(fes1,fes2);
% Shift downwards to center the geometry
fens = transform_apply(fens,@(x,d)x+[0,0,-Length/2], []);

% Convert to triangles
[fens,fes] = Q4_to_T3(fens,fes,struct( 'other_dimension', 1 ));
% Retrieve output arrays
xyz=fens.xyz;
conn=fes.conn;
% Produce a plot
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');
view (3)
bdfes=mesh_boundary(fes,struct('other_dimension', 100 ));
if (count (bdfes)>0)
    error('Inconsistent mesh: merging failed');
end
%  drawmesh({fens,bdfes},'fes', 'facecolor', 'y')
end