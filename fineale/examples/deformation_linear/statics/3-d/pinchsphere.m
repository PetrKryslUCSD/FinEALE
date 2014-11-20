% The pinched hemisphere benchmark
% The reference below states:
%
% The spherical shell shown in Fig. 9 is our proposed doubly-curved shell
% problem. Note that the equator is a free edge so that the problem
% represents a hemisphere with four point loads alternating in sign at 90 °
% intervals on the equator. The hole at the top has been introduced to
% avoid the use of triangles ne:~r the axis of revolution. Convergence can
% be studied by varying  mesh size. Both membrane and bending strains
% contribute significantly to the radial displacement at the load point. A
% theoretical value of the displacement under load has been computed for a
% slightly different configuration [7] in which the hole at the axis is
% closed.
%
%     In our example the shell is closed at the top, and we avoid the use
%     of triangles on the axis of symmetry..
%
% Macneal RH, Harder RL (1985) A proposed standard set of problems to test
% finite element accuracy. Finite Elements in Analysis and Design 1: 3-20.
function gv =pinchsphere
E=6.825e7;
nu=0.3;
thickness = 0.04;
% analytical solution for the vertical deflection under the load
analyt_sol=0.0924;
graphics = ~true;
R=10;
u_scale = 20;
nlayers = 1;
ns=1:5;

clc;
disp(['nlayers = ' num2str(nlayers) ';']);
tolerance  =thickness/nlayers/1000;

% Material
prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;


eltyd(eix).description ='H20R';
eltyd(eix).mf =@H8_to_H20;
eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;
%
eltyd(eix).description ='H20';
eltyd(eix).mf =@H8_to_H20;
eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

eltyd(eix).description ='H27';
eltyd(eix).mf =@H8_to_H27;
eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

eltyd(eix).description ='H20-Bbar';
eltyd(eix).mf =@H8_to_H20;
eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
    'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
    'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
    'pv_bfun',@(p)[1;p(1);p(2);p(3)],'nconstrained',1));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;


legends ={};
for eix = 1:length(eltyd)
    
    uzs =[];
    for n  = ns
        
        % Mesh
        [fens,fes]=H_spherical_shell_mesh(R-thickness/2,thickness,n,nlayers, eltyd(eix).mf);
        %         count(fens)
        %             drawmesh({fens,fes},'fes','facecolor','red')
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.fes= fes;
        region.femm= eltyd(eix).femmf(fes);
        model_data.region{1} =region;
        
        clear essential
        essential.component= [2];
        essential.fixed_value= 0;
        essential.node_list = fenode_select (fens,struct ('box',[-10000 10000 0 0 -10000 10000],'inflate',tolerance));
        model_data.boundary_conditions.essential{1} = essential;
        
        clear essential
        essential.component= [1];
        essential.fixed_value= 0;
        essential.node_list = fenode_select (fens,struct ('box',[0 0 -10000 10000 -10000 10000],'inflate',tolerance));
        model_data.boundary_conditions.essential{2} = essential;
        
        clear essential
        essential.component= [];
        essential.fixed_value= 0;
        essential.node_list = fenode_select (fens,struct ('box',[0 0 0 0 R-thickness/2 R-thickness/2],'inflate',tolerance));
        model_data.boundary_conditions.essential{3} = essential;
        
        clear nodal_force
        nodal_force.node_list =fenode_select (fens,struct ('box',[-inf inf 0 0 0  0],'inflate',tolerance));
        nodal_force.force= [1/length(nodal_force.node_list);0;0];
        model_data.boundary_conditions.nodal_force{1} = nodal_force;
        
        clear nodal_force
        nodal_force.node_list =fenode_select (fens,struct ('box',[0 0 -inf inf 0  0],'inflate',tolerance));
        nodal_force.force= [0;-1/length(nodal_force.node_list);0];
        model_data.boundary_conditions.nodal_force{2} = nodal_force;
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        corn=fenode_select (fens,struct('box',[0 0 -inf inf 0  0],'inflate',tolerance));
        ucorn= gather_values(model_data.u, corn);
        uz=mean(ucorn(:,2));
        nd=abs(uz)/analyt_sol;
        disp(['%' eltyd(eix).description ': ' 'Deflection under the load: ' num2str(nd*100) '%'])
        uzs =[uzs , uz];
        
        if (graphics )
            clear options
            options.u_scale= u_scale;
            options=deformation_plot_deformation(model_data, options);
        else
            
        end
    end
    if (graphics )
    else
        %             semilogx(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
        plot(2*2.^ns,abs(uzs/analyt_sol),name_to_style(eltyd(eix).description),'linewidth',3);  hold on
        title(['Pinched spherical shell, nlayers  =' num2str(nlayers)])
        figure (gcf); pause (1)
    end
    legends{end+1} =eltyd(eix).description;
    format long
    disp(['Data{end+1}=[']);
    disp([2*2.^ns;uzs]);
    disp(['];' 'Style{end+1}='''  name_to_style(eltyd(eix).description) ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
end
legend (legends);
labels(  'Number of equations', 'Estimated true error')
grid on
set_graphics_defaults
end