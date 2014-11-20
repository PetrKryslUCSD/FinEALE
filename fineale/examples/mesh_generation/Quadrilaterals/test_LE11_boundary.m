% Extract the boundary of the mesh of the LE11 benchmark.
function test_LE11_boundary
  pu= physical_units_struct;
  X=[1.    , 0.;%A
    1.4   , 0.;%B
    0.995184726672197   0.098017140329561;
    1.393258617341076 0.137223996461385;
    0.980785,0.195090;%
    1.37309939,0.27312645;
    0.956940335732209   0.290284677254462
    1.339716470025092 0.406398548156247
    0.9238795, 0.38268;%C
    1.2124, 0.7;%D
    0.7071, 0.7071;%E
    1.1062, 1.045;%F
    0.7071, (0.7071+1.79)/2;%(E+H)/2
    1.    , 1.39;%G
    0.7071, 1.79;%H
    1.    , 1.79;%I
    ]*pu.M;
  fens=fenode_set(struct('xyz',X))
  fes=fe_set_Q4(struct('conn',[1,2,4,3;3,4,6,5;5,6,8,7;7,8,10,9;9,10,12,11;11,12,14,13;13,14,16,15]))
  % [fens,fes]=Q4_refine(fens,fes);
  length_units=pu.MM;
  options.thickness =1.0;
  gv=drawmesh({fens,fes},'fes','facecolor',[1,1,1]/1.2,'length_units',length_units);
  bfes=mesh_boundary(fes,[]);
  gv=drawmesh({fens,bfes},'gv',gv,'fes','facecolor','none','edgecolor','m','linewidth',2,'length_units',length_units);
  bfes=mesh_boundary(subset(bfes,3),[]);
  gv=drawmesh({fens,bfes},'gv',gv,'fes','facecolor','r','linewidth',2,'length_units',length_units);
  draw_arrow(gv,X(1,:),X(13,:)-X(1,:),struct('color','y','length_units',length_units));
  draw_axes(gv, struct('length',0.5,'length_units',length_units))
  draw_text(gv, [0.6,0.4], 'Here', struct('color','c','length_units',length_units))
  draw_ellipsoid(gv, [0,1], eye(3), [1,1,1]./[2,3,4], struct('tessel',16,'color','c','length_units',length_units))
draw_cylinder(gv, X(end,:), X(end-1,:), 0.1, 0.08, struct('tessel',16,'facecolor','g','length_units',length_units))
 draw_planes(gv, struct('length',0.3,'facealpha',0.3,'length_units',length_units))
  headlight(gv)
  view(2)
 