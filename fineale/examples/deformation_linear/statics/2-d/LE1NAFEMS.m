% LE1 NAFEMS benchmark
function f
pu= physical_units_struct;
% Parameters:
E = 210e3*pu.MEGA*pu.PA;% 210e3 MPa
nu = 0.3;
p = 10*pu.MEGA*pu.PA;% 10 MPA Outward pressure on the outside ellipse
sigma_yD= 92.7*pu.MEGA*pu.PA;% tensile stress at [2.0, 0.0] meters
n=12;% number of elements through the thickness

integration_rule =gauss_rule(struct('dim', 2,'order', 2));
Edge_integration_rule =gauss_rule(struct('dim', 1,'order', 2));
Convertf=[];
Style ='ks-'; Label='Q4';

% integration_rule = gauss_rule(struct('dim', 2,'order', 3));
% Edge_integration_rule =gauss_rule(struct('dim', 1,'order', 3));
% Convertf=@Q4_to_Q16;
% Style ='ks-'; Label='Q16';
%
% integration_rule = gauss_rule(struct('dim', 2,'order', 2));
% Edge_integration_rule =gauss_rule(struct('dim', 1,'order', 3));
% Convertf=@Q4_to_Q8;
% Style ='ks-'; Label='Q8';
%
% integration_rule = simpson_1_3_rule(struct('dim',2));
% Edge_integration_rule =simpson_1_3_rule(struct('dim',1));
% Convertf=@Q4_to_Q8;
% Style ='ks-'; Label='Q8';

% 
% Mesh'
[fens,fes]=Q4_block(1.0,pi/2, n, n*2, struct('other_dimension', 0.1*1000));
if (~isempty( Convertf ))
  [fens,fes] = Convertf(fens,fes, struct('other_dimension', 0.1*1000));
end
bdry_fes = mesh_boundary(fes, struct('other_dimension', 0.1*1000));
icl = fe_select(fens, bdry_fes, struct('box', [1,1,0,pi/2],'inflate',1/100));
xy=fens.xyz;
for i=1:count(fens)
    t=xy(i,1); a=xy(i,2);
    xy(i,:)=1000*[(t*3.25+(1-t)*2)*cos(a) (t*2.75+(1-t)*1)*sin(a)];
end
fens.xyz=xy;




    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.E =E;
    region.nu =nu;
    region.reduction ='stress';
    region.fes= fes;
    region.integration_rule = integration_rule;
    region.Rm =[];
    model_data.region{1} =region;
    
    clear essential
    essential.component= 2;
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0 Inf 0 0],'inflate',1/10000));
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.component= 1;
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0 0 0 Inf],'inflate',1/10000));
    model_data.boundary_conditions.essential{2} = essential;
    
    clear traction
    traction.fes= subset(bdry_fes,icl);
    traction.integration_rule = Edge_integration_rule;
    traction.traction = @(x) (p*[2.75/3.25*x(1),3.25/2.75*x(2)]'/norm([2.75/3.25*x(1),3.25/2.75*x(2)]));
    model_data.boundary_conditions.traction{1} = traction;
    
    % Solve
    model_data =deformation_linear_statics(model_data);
    
% toc, tic
% Plot
gv=graphic_viewer;
gv=reset (gv,struct ('limits',1000*[0 1.06*3.25 0 1.6*3.25]));
scale=1000;
cmap = jet;
fld = field_from_integration_points(model_data.region{1}.femm, model_data.geom, model_data.u, 0*model_data.u, 0, [], 'Cauchy',2);
corner=fenode_select (fens,struct('box',1000*[2.0 2.0 0 0],'inflate',1/10000));
disp( ['Corner displacement=' num2str(gather_values (model_data.u, corner)) ' M'])
disp( ['Corner stress sigma_y=' num2str(gather_values (fld, corner)/(pu.MEGA*pu.PA)) ' MPa'])
disp( [' i.e. ' num2str(gather_values (fld, corner)/sigma_yD*100) '% of the Reference value'])
nvals=fld.values;%min(nvals),max(nvals)
nvalsrange=[min(nvals),max(nvals)];
dcm=data_colormap(struct ('range',nvalsrange, 'colormap',cmap));
colorfield=nodal_field(struct ('name', ['colorfield'], 'data',map_data(dcm, nvals)));
draw(model_data.region{1}.femm,gv, struct ('x', model_data.geom,  'u', +scale*model_data.u,'colorfield',colorfield, 'shrink',1.0));
draw(model_data.region{1}.femm,gv, struct ('x', model_data.geom,  'u', +0*model_data.u,'facecolor','none', 'shrink',1.0));
% draw(efemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','red', 'shrink',1.0));
colormap(cmap);
draw_colorbar(gv,struct('position',[0.75 0.15 0.025 0.7],'minmax',nvalsrange,'label','\sigma_y'));
view (2)
set_graphics_defaults
% assignin('caller','fineale_test_passed',((norm([93.029698724563474]-gather(fld, corner,'values')))<1e-9))