% Slot in a rubber block under plane strain.
% @article{
%    author = {Nakshatrala, K. B. and Masud, A. and Hjelmstad, K. D.},
%    title = {On finite element formulations for nearly incompressible linear elasticity},
%    journal = {Computational Mechanics},
%    volume = {41},
%    number = {4},
%    pages = {547-561},
%    ISSN = {0178-7675},
%    DOI = {10.1007/s00466-007-0212-8},
%    year = {2008}
% }

function  Slot_T10
% Parameters:
nu =0.49995;
E=2.4e6;
tol=1/10000;
graphics = true;
u_scale=1;
%     nass={'Slot-coarse.nas','Slot-fine.nas','Slot-extra-fine.nas'};
nass={'Slot-coarser.nas'};

ns=1;

format long
clc
disp(['% ' mfilename]);
disp(['nu=' num2str(nu,12)]);

prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;

%     eltyd(eix).description ='T10';% tetrahedron
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eltyd(eix).styl='k^--';
%     eix=eix+1;
%
    eltyd(eix).description ='T10-SRI';
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
        'integration_rule_volumetric',tet_rule(struct('npts',1)),...
        'integration_rule_deviatoric',tet_rule(struct('npts',4))));
    eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    eltyd(eix).styl='kv-';
    eix=eix+1;


for eix = 1:length(eltyd)
    uzs=[];
    nfreedofs = [];
    
    for nas=nass  %,16,24,32,48
        % Create the mesh and initialize the geometry
        [fens,fes]= nastran_mesh_import(nas{1});
        [fens,fes] = T10_refine(fens,fes);
        %             drawmesh({fens,fes},'fes'); labels
        
        bfes=mesh_boundary(fes,[]);
        z1=fe_select(fens,bfes,struct('facing',true, 'direction',[0,0,-1],'tolerance',0.999));
        z2=fe_select(fens,bfes,struct('facing',true, 'direction',[0,0,+1],'tolerance',0.999));
        x0=fe_select(fens,bfes,struct('box',[0,0,-inf,inf,-inf,inf],'inflate',tol));
        x1=fe_select(fens,bfes,struct('box',[1,1,-inf,inf,-inf,inf],'inflate',tol));
        y0=fe_select(fens,bfes,struct('box',[-inf,inf,0,0,-inf,inf],'inflate',tol));
        y1=fe_select(fens,bfes,struct('box',[-inf,inf,1,1,-inf,inf],'inflate',tol));
        outlineprop=struct( 'other_dimension' ,0.5);
        outline=mesh_boundary(subset(bfes, [z1]),outlineprop);
        outline=cat(outline,mesh_boundary(subset(bfes, [z2]),outlineprop));
        outline=cat(outline,mesh_boundary(subset(bfes, [x0]),outlineprop));
        outline=cat(outline,mesh_boundary(subset(bfes, [x1]),outlineprop));
        outline=cat(outline,mesh_boundary(subset(bfes, [y0]),outlineprop));
        outline=cat(outline,mesh_boundary(subset(bfes, [y1]),outlineprop));
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.fes= fes;
        region.femm= eltyd(eix).femmf(fes);
        model_data.region{1} =region;
        
        clear essential % Symmetry plane and the constrained side surface
        essential.component= [1];
        essential.fixed_value= 0;
        essential.node_list = [fenode_select(fens,struct ('box',[0,0,-inf,inf,-inf,inf],'inflate',tol)),...
            fenode_select(fens,struct ('box',[1,1,-inf,inf,-inf,inf],'inflate',tol))];
        model_data.boundary_conditions.essential{1} = essential;
        
        clear essential % Symmetry plane
        essential.component= [2];
        essential.fixed_value= 0;
        essential.node_list = fenode_select (fens,struct ('box',[-inf,inf,0,0,-inf,inf],'inflate',tol));
        model_data.boundary_conditions.essential{2} = essential;
        
        clear essential % Enforcement of plain strain conditions
        essential.component= [3];
        essential.fixed_value= 0;
        essential.node_list = [fenode_select(fens,struct('box',[-inf,inf,-inf,inf,-0.1/2,-0.1/2],'inflate',tol)),...
            fenode_select(fens,struct('box',[-inf,inf,-inf,inf,0.1/2,0.1/2],'inflate',tol))];
        model_data.boundary_conditions.essential{3} = essential;
        
        clear essential % Displaced support section
        essential.component= [2];
        essential.fixed_value= 0.01;
        essential.node_list = fenode_select (fens,struct ('box',[-inf,inf,1,1,-inf,inf],'inflate',tol));
        model_data.boundary_conditions.essential{4} = essential;
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        if (graphics )
            %                 model_data.postprocessing.u_scale= u_scale;
            %                 model_data=deformation_plot_deformation(model_data);
            figure
            model_data.postprocessing.u_scale= u_scale;
            model_data.postprocessing.output='pressure';
            model_data.postprocessing.camera=[ -3.8327    2.5419    7.2914    0.4406    0.5615   -0.0758    0.1153    0.9740   -0.1950    6.5627];
            model_data.postprocessing.stress_component=1;
            model_data.postprocessing.boundary_only=true;
            model_data.postprocessing.stress_range=[0,800000];
            model_data.postprocessing.map_to_color_fun =@map_to_color_fun;
            model_data.postprocessing.add_decorations =@add_decorations;
            m=parula; %m=m(end:-1:1,:);
            model_data.postprocessing.cmap=m;
            model_data.postprocessing.add_to_scene=@add_to_scene;
            model_data.postprocessing.edgecolor='none';
            model_data=deformation_plot_stress_elementwise(model_data);
            %                 model_data=deformation_plot_stress(model_data);
            %               headlight(model_data.postprocessing.gv);
        end
        
    end
    
end

    function gv=add_to_scene(gv,~)
        %         title(eltyd(eix).description);
        draw(outline, gv, struct ('x',model_data.geom, 'u',u_scale*model_data.u,...
            'edgecolor','k'));
    end

minmin=inf;
maxmax=-inf;
    function [colorfield,dcm]=map_to_color_fun(fld,dcm);
        fld=(-1)*fld;
        minmin=min(fld.values);
        maxmax=max(fld.values);
        %         dcm=data_colormap(struct('range',[minmin,maxmax],'colormap',parula));
        colorfield=nodal_field(struct ('name', ['colorfield'], 'data',...
            map_data(dcm, fld.values)));
    end

    function gv=add_decorations(gv,dcm);
        %         xlabel('X')
        %         ylabel('Y')
        %         zlabel('Z')
        axis off
        draw_colorbar(gv,struct('colormap',dcm.colormap,...
            'position',[0.81, 0.1, 0.025, 0.5],...
            'minmax', dcm.range,...
            'label',['Pressure $p$'], 'fontname', 'Times', 'interpreter', 'latex'));
        
    end
end