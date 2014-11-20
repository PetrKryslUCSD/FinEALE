% Generation of a spherical cavity mesh.
R=5.0;% radius of the sphere, millimeters
nR=3;
Ro=nR*R; % Outer radius, millimeters
nperR=6;%
nlayers=nperR;

tol=R/nperR/1000;

% % %     % Mesh

[fens,fes]=H8_sphere_n(R,nperR);

bdry_fes=mesh_boundary(fes);
clear options
options. facing=true;
options. direction=[1,1,1];
Intl=fe_select(fens, bdry_fes, options);
clear fes
[fens,fes] = H8_extrude_Q4(fens,subset(bdry_fes,Intl),nlayers,@(xyz, layer)(R+layer/nlayers*(Ro-R))*xyz/norm(xyz));
connected = find_unconn_fens(fens, fes);
[fens, new_numbering] =compact_fens(fens, connected);
fes = update_conn(fes, new_numbering);


bdry_fes=mesh_boundary(fes);

gv=drawmesh({fens,bdry_fes},'fes','facecolor','red');
