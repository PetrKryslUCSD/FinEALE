% MEMS actuator.   Thermal analysis.
m=physical_units_machine;
x0= 0;
x1=x0+5e-6;
x2=x1+10e-6;
x3=x2+10e-6;
x4=x3+10e-6;
y0=0;
y4=250e-6;
y3=y4-10e-6;
y2=y3-10e-6;
y1=y2-10e-6;
t=2e-6;
h=0.1e-6;
z0=0;
z3=2*t+h;
z2=z3-t;
z1=z2-h;
m1=1*2;
m2=1*2;
m3=1*2;
m4=1*2;
n1=2*2;
n2=1*2;
n3=1*2;
n4=2*2;
n5=2*2;
p1=1*2;
p2=1;
p3=1*2;
kappa=157*eye(3, 3); % W/m/K, conductivity matrix
DV=5;% voltage drop in volt
l =2*(y1+y2)/2+2*(x1+x2)/2;% length of the conductor
resistivity = 1.1e-5;% Ohm m
Q=DV^2/resistivity/l^2;% rate of Joule heating, W/m^3
T_substrate=293;% substrate temperature in degrees Kelvin

[fens,fes] = H8_hexahedron([x1,y0,z0;x2,y1,z1],m2,n1,p1);
[fens1,fes1] = H8_hexahedron([x1,y1,z0;x2,y2,z1],m2,n2,p1);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x0,y1,z0;x1,y2,z1],m1,n2,p1);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x0,y1,z1;x1,y2,z2],m1,n2,p2);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x0,y1,z2;x1,y2,z3],m1,n2,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x0,y2,z2;x1,y3,z3],m1,n3,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x0,y3,z2;x1,y4,z3],m1,n4,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x1,y3,z2;x3,y4,z3],m4,n4,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x3,y3,z2;x4,y4,z3],m3,n4,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
[fens1,fes1] = H8_hexahedron([x3,y0,z2;x4,y3,z3],m3,n5,p3);
[fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, eps);
fes= cat(fes1,fes2);
% gv=drawmesh({fens,fes},'fes','facecolor','red');
[fens,fes] = H8_to_H20(fens,fes);
count(fens)

hotprop=property_heat_diffusion(struct('thermal_conductivity',kappa,'source',Q));
hotmater=material_heat_diffusion (struct('property',hotprop));
coldprop=property_heat_diffusion(struct('thermal_conductivity',kappa,'source',0.0));
coldmater=material_heat_diffusion (struct('property',coldprop));
cl= fe_select(fens, fes, struct('box',[x0,x2,y0,y2,z0,z1],'inflate',t/100));
hotfemm = femm_heat_diffusion (struct ('material',hotmater,...
    'fes',subset(fes,cl),...
    'integration_rule',gauss_rule(struct( 'dim',3,'order',2))));
coldfemm = femm_heat_diffusion (struct ('material',coldmater,...
    'fes',subset(fes,setdiff((1:count(fes)),cl)),...
    'integration_rule',gauss_rule(struct( 'dim',3,'order',2))));
geom = nodal_field(struct('name',['geom'], 'dim', 3, 'fens',fens));
theta=nodal_field(struct('name',['theta'], 'dim', 1, 'nfens',geom.nfens));
fenids=fenode_select(fens,struct('box',[x0,x4,y0,y0,z0,z3],...
    'inflate', t/1000)) ; % fixed temperature on substrate
fixed=fenids*0+1;
comp=[];
val=zeros(length(fenids),1)+T_substrate;
theta = set_ebc(theta, fenids, fixed, comp, val);
theta = apply_ebc (theta);

theta = numberdofs (theta);
K = conductivity(hotfemm, sysmat_assembler_sparse, geom, theta)+...
    conductivity(coldfemm, sysmat_assembler_sparse, geom, theta);
fi= force_intensity(struct('magn',Q));
F = distrib_loads(hotfemm, sysvec_assembler, geom, theta, fi, 3)+...
 nz_ebc_loads_conductivity(hotfemm, sysvec_assembler, geom, theta)+...
    nz_ebc_loads_conductivity(coldfemm, sysvec_assembler, geom, theta);
theta = scatter_sysvec(theta, K\F);

y_i=gather_values(geom,fenode_select(fens,struct('box',[x1,x1,y0,y1,z1,z1],...
    'inflate', t/1000)));
T_i=gather_values(theta,fenode_select(fens,struct('box',[x1,x1,y0,y1,z1,z1],...
    'inflate', t/1000)));
[ignore,ix]=sort(y_i(:,2));
figure;
plot(y_i(ix,2),T_i(ix),'bd','linewidth',3); hold on
y_o=gather_values(geom,fenode_select(fens,struct('box',[x3,x3,y0,y3,z2,z2],...
    'inflate', t/1000)));
T_o=gather_values(theta,fenode_select(fens,struct('box',[x3,x3,y0,y3,z2,z2],...
    'inflate', t/1000)));
[ignore,ix]=sort(y_o(:,2));
plot(y_o(ix,2),T_o(ix),'bo','linewidth',3); hold off

% Plot
gv=graphic_viewer;
gv=reset (gv,[]);
T=theta.values;
cmap= hot;
dcm=data_colormap(struct('range',[min(T),max(T)],'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
    map_data(dcm, T)));
draw(mesh_boundary (hotfemm.fes, []), gv, struct ('x',geom, 'u',0*geom,...
    'colorfield',colorfield, 'shrink',1));

draw(mesh_boundary (coldfemm.fes, []), gv, struct ('x',geom, 'u',0*geom,...
    'colorfield',colorfield, 'shrink',1));
draw_colorbar(gv, struct('colormap',cmap,'label','Temperature',...
    'minmax',[min(T),max(T)]))
interact(gv);


export_to_abaqus= true; SurfElType ='SFM3D4';
        ElType ='C3D8';
        ElType ='C3D8H';
                        ElType ='C3D8I';
        %             ElType ='C3D8IH';
        %                         ElType ='C3D8RH';
        ElType ='C3D20R'; SurfElType ='SFM3D8';;
    
    if (export_to_abaqus)
                date_now = clock;
                s = strcat(num2str(date_now(1)),'-',num2str(date_now(2)),'-', num2str(date_now(3)), '-',num2str(date_now(4)), '-',num2str(date_now(5)));
                AE=Abaqus_exporter;
                AE.open([mfilename '-' ElType '-'  s '.inp']);
                AE.HEADING([mfilename ' ' 'ElType=' ElType ]);
                AE.PART('part1');
                AE.END_PART();
                AE.ASSEMBLY('ASSEM1');
                AE.INSTANCE('INSTNC1','PART1');
                AE.NODE(1000*geom.values);
                AE.ELEMENT(ElType,'All',1,fes.conn);
                AE.END_INSTANCE();
                AE.END_ASSEMBLY();
                AE.close();
            end