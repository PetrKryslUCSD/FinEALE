% Generate the surface of a submarine.
function test_submarine_surface1
    u=physical_units_struct;
% The axis of the submarine hull is the X axis
% The tower axis is the Z axis
% radius of the  hull
R_h=5000*u.MM;
nperR_h=4*2;
% major radius of the tower
R_t=0.8*R_h;
nperR_t=nperR_h;
% minor radius of the tower
r_t=0.25*R_h;
nperr_t=nperR_t;
% tower height
L_t=0.5*R_h;
nperL_t=4;
% Fore length
L_f=2.5*R_h;
nperL_f=nperR_h;
% Aft length
L_a=4.5*R_h;
nperL_a=nperR_h;
% Cone length
L_c=3.5*R_h;
nperL_c=nperR_h;
% Cone radius
r_c=.3*R_h;
% Top fin
L_t=.5*R_h;
L_F=0.8*R_h;
nperL_F=2;
H_T=1.25*R_h;
nperH_T=2;
H_B=1.25*R_h;
nperH_B=2;

% geometrical tolerance
tol=r_t/nperr_t/2;

% Initialize the graphic model
gv=graphic_viewer;
 gv=reset(gv,[]);
 
% fore section
[fensf,fef]=Q4_elliphole(R_t,r_t,L_f+R_t,pi*R_h,nperL_f,nperR_h,nperR_h,[]);
% aft section
[fensa,fea]=Q4_elliphole(R_t,r_t,L_a+R_t,pi*R_h,nperL_a,nperR_h,nperR_h,[]);
[fensa,fea] = mirror_mesh(fensa, fea, [-1,0], [0,0], @(c)c([1, 4, 3, 2]));
% Merge fore and aft sections
[fens,fe1,fe2] = merge_meshes(fensf, fef, fensa, fea, tol);
fe=cat(fe1,fe2);
% Expand the coordinates to 3-D
xy=fens.xyz; fens.xyz=[xy,zeros(size(xy,1),1) ];
% Roll up the hull
fens = transform_apply(fens,@(x,d)[x(1),R_h*sin(x(2)/R_h),R_h*cos(x(2)/R_h)], []);
% Mirror the hull
[fens1,fe1] = mirror_mesh(fens, fe, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fe1,fe2] = merge_meshes(fens1, fe1, fens, fe, tol);
fe=cat(fe1,fe2);
% Save the hull
fensh =fens;
feh=fe;
% Render the hull
% gv=drawmesh({fensh,feh},'gv',gv,'fe', 'facecolor', 'y');
 

% Make 1/8 of the spherical surface
[fens,fe]=Q4_sphere_n(R_h,nperR_h/2,1);
%     Mirror to get 1/4
[fens1,fe1] = mirror_mesh(fens, fe, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fe1,fe2] = merge_meshes(fens1, fe1, fens, fe, tol);
fe=cat(fe1,fe2);
%     Mirror again to get 1/2: this is the nose
[fens1,fe1] = mirror_mesh(fens, fe, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fe1,fe2] = merge_meshes(fens1, fe1, fens, fe, tol);
fe=cat(fe1,fe2);
% Rotate the nose
fens = transform_apply(fens,@(x,d)[x(3),x(2),-x(1)], []);
% Translate the nose
fens = transform_apply(fens,@(x,d)x+[R_t+L_f,0,0], []);
% Save the nose
fensn =fens;
fen=fe;
% gv=drawmesh({fensn,fen},'gv',gv,'fe', 'facecolor', 'r');
% Make the cone
[fens,fe] = Q4_block(L_c,pi*R_h,nperL_c,nperR_h,struct( 'other_dimension', 1 ));
% Expand the coordinates to 3-D
xy=fens.xyz; fens.xyz= [xy,zeros(size(xy,1),1) ];
% Roll up the cone
fens = transform_apply(fens,@(x,d)[x(1),R_h*sin(x(2)/R_h),R_h*cos(x(2)/R_h)], []);
% taper the cone
fens = transform_apply(fens,@(x,d)[x(1),(1-x(1)/L_c*(1-R_h/r_c))*r_c/R_h*[x(2),x(3)]], []);
% Translate the cone
fens = transform_apply(fens,@(x,d)x+[-L_c-L_a-R_t,0,0], []);
% Close the cone
[fensC,feC]=Q4_circle_n(r_c,nperR_h/2,1.0);
% Expand the coordinates to 3-D
xy=fensC.xyz; fensC.xyz= [xy,zeros(size(xy,1),1) ];
%     Mirror to get 1/2
[fens1,fe1] = mirror_mesh(fensC, feC, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fensC,fe1,fe2] = merge_meshes(fens1, fe1, fensC, feC, tol);
feC=cat(fe1,fe2);
% Rotate the closure of the cone
fensC = transform_apply(fensC,@(x,d)[x(3),x(2),-x(1)], []);
% Translate the closure of the cone
fensC = transform_apply(fensC,@(x,d)x+[-L_c-L_a-R_t,0,0], []);
% merge cone with closure
[fens,fe1,fe2] = merge_meshes(fensC, feC, fens, fe, tol);
fe=cat(fe1,fe2);
%     Mirror to get the whole
[fens1,fe1] = mirror_mesh(fens, fe, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fens,fe1,fe2] = merge_meshes(fens1, fe1, fens, fe, tol);
fe=cat(fe1,fe2);
% Save the cone
fensc =fens;
fec=fe;
% gv=drawmesh({fensc,fec},'gv',gv,'fe', 'facecolor', 'c');

% Merge everything together
fens=fensh; fe=feh;
tol=R_h/nperR_h/2;
[fens,fe1,fe2] = merge_meshes(fensn, fen, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensc, fec, fens, fe, tol);
fe=cat(fe1,fe2);

% Make the tower
tbdfe=mesh_boundary(fe,struct('other_dimension', 1));
% gv=drawmesh({fens,tbdfe},'gv',gv,'fe', 'facecolor', 'r');
fenst=fens;
[fens,fe] = Q4_extrude_L2(fenst,tbdfe,nperL_t,@(x, layer)x+[0,0,layer/nperL_t*L_t]);
% Save the Tower
fenst =fens;
fet=fe;
% gv=drawmesh({fenst,fet},'gv',gv,'fe', 'facecolor', 'm');

% Close the tower
[fensT,feT]=Q4_circle_n(r_t,2*nperR_h,1.0);
% Expand the coordinates to 3-D
xy=fensT.xyz; fensT.xyz= [xy,zeros(size(xy,1),1) ];
% Mirror
tol=r_t/nperR_h/5;
[fens1,fe1] = mirror_mesh(fensT, feT, [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fensT,fe1,fe2] = merge_meshes(fens1, fe1, fensT, feT, tol);
feT=cat(fe1,fe2);
[fens1,fe1] = mirror_mesh(fensT, feT, [0,-1,0], [0,0,0], @(c)c([1, 4, 3, 2]));
[fensT,fe1,fe2] = merge_meshes(fens1, fe1, fensT, feT, tol);
feT=cat(fe1,fe2);
% Make into an ellipse
fensT = transform_apply(fensT,@(x,d)[(R_t/r_t)*x(1),x(2),x(3)], []);
% Roll up 
fensT = transform_apply(fensT,@(x,d)[x(1),R_h*sin(x(2)/R_h),R_h*cos(x(2)/R_h)], []);
% Translate
fensT = transform_apply(fensT,@(x,d)x+[0,0,L_t], []);
% gv=drawmesh({fensT,feT},'gv',gv,'fe', 'facecolor', 'm');


% Make the fin
% Top part
x1= [-(R_t+L_a+L_c)+L_t,0,r_c+(R_h-r_c)/L_c*L_t];
x2= x1+[0,0,H_T-x1(3)];
x3=x2+[L_F,0,0];
x4= [-(R_t+L_a+L_c)+L_t+L_F,0,r_c+(R_h-r_c)/L_c*(L_t+L_F)];
[fens,fe] = Q4_quadrilateral([x1;x2;x3;x4],nperL_F,nperH_T,struct('other_dimension', 1));
% Save the Top Fin
fensFT =fens;
feFT=fe;
% gv=drawmesh({fensFT,feFT},'gv',gv,'fe', 'facecolor', 'b');
% Bottom part
x1= [-(R_t+L_a+L_c)+L_t,0,-(r_c+(R_h-r_c)/L_c*L_t)];
x2= x1+[0,0,-H_B-x1(3)];
x3=x2+[L_F,0,0];
x4= [-(R_t+L_a+L_c)+L_t+L_F,0,-(r_c+(R_h-r_c)/L_c*(L_t+L_F))];
[fens,fe] = Q4_quadrilateral([x1;x2;x3;x4],nperL_F,nperH_B,struct('other_dimension', 1));
% Save the Bottom Fin
fensFB =fens;
feFB=fe;
% gv=drawmesh({fensFB,feFB},'gv',gv,'fe', 'facecolor', 'b');

% merge everything, but don't glue together the vertices
tol=0;
fens=fensh; fe=feh;
[fens,fe1,fe2] = merge_meshes(fensn, fen, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensc, fec, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensc, fec, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fenst, fet, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensFT, feFT, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensFB, feFB, fens, fe, tol);
fe=cat(fe1,fe2);
[fens,fe1,fe2] = merge_meshes(fensT, feT, fens, fe, tol);
fe=cat(fe1,fe2);

% Convert to triangles
% [fens,fe] = Q4_to_T3_sd(fens,fe,struct( 'other_dimension', 1 ));
% Retrieve output arrays
xyz=fens.xyz;
conn=fe.conn;
% Produce a plot
h=patch('Faces', conn, ...
    'Vertices', xyz,'BackFaceLighting','lit',...
    'AmbientStrength',0.75,...
    'FaceColor', 'r','EdgeColor','k');


 labels('X','Y','Z');
 interact(gv,struct('snap_points',xyz));
end