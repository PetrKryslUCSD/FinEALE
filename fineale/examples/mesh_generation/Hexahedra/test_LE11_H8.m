% Test generation of hexahedral meshfor the LE11 benchmark
function test_LE11
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
 nLayers=1;
 angslice =pi/16;
[fens,fes] = H8_extrude_Q4(fens,fes,nLayers,@(xyz,k)[xyz,0]+(k-1/2)/nLayers*[0,0,angslice]);
bfes=mesh_boundary(fes);
 f1l=fe_select (fens,bfes,struct('box',[-inf,inf,-inf,inf,-angslice/2,-angslice/2],'inflate',100*eps));
 f2l=fe_select (fens,bfes,struct('box',[-inf,inf,-inf,inf,angslice/2,angslice/2],'inflate',100*eps));
fens = transform_apply(fens,@(x,d)[x(1)*cos(x(3)),x(2),-x(1)*sin(x(3))], []);

    gv=drawmesh({fens,fes},'fes','facecolor','none');
 gv=drawmesh({fens,subset(bfes,f1l)},'gv',gv,'fes','facecolor','y');
gv=drawmesh({fens,subset(bfes,f2l)},'gv',gv,'fes','facecolor','r');

 
 
 