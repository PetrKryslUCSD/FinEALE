% The pinched hemisphere benchmark, tetrahedral elements
function gv =pinchsphere_T10
    disp('Pinched sphere');
    E=6.825e7;
    nu=0.3;
    thickness = 0.04;
    tolerance  =thickness/1000;
    % analytical solution for the vertical deflection under the load
    analyt_sol=0.0924;
    graphics = ~true;
    R=10;
    u_scale = 20;
    
    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    meshes =  {'hemisphere000.inp','hemisphere00.inp','hemisphere0.inp'};
    clear eltyd
    eix=1;
    
%     eltyd(eix).description ='NICE-GSRI-T10';% tetrahedron
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_nice_gsri(struct('fes',fes,'material',mater,...
%         'integration_rule_constrained',simplex_nq_rule(struct('dim',3, 'order',2)),...
%         'integration_rule_unconstrained',tet_rule(struct('npts',4))));
%     eltyd(eix).styl='c^-';
%     eix=eix+1;
    
        eltyd(eix).description ='T10';% tetrahedron
        eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
            'integration_rule',tet_rule(struct('npts',4))));
        eltyd(eix).styl='gd-';
        eix=eix+1;
    
    %     eltyd(eix).description ='T10-nq';% tetrahedron
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
    %         'integration_rule',simplex_nq_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).styl='rd--';
    %     eix=eix+1;
    
%         eltyd(eix).description ='T10-SRI';
%         eltyd(eix).mf =@T10_block;
%         eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%             'integration_rule_volumetric',tet_rule(struct('npts',1)),...
%             'integration_rule_deviatoric',tet_rule(struct('npts',5)),'split', 'bulk'));
%         eltyd(eix).styl='gp-';
%         eix=eix+1;
%     
        eltyd(eix).description ='T10-GSRI';
        eltyd(eix).mf =@T10_block;
        eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
            'integration_rule_constrained',tet_rule(struct('npts',1)),...
            'integration_rule_unconstrained',tet_rule(struct('npts',4))));
        eltyd(eix).styl='r*-';
        eix=eix+1;
    
    %     eltyd(eix).description ='T10-Bbar';
    %     eltyd(eix).mf =@T10_block;
    %     eltyd(eix).femmf =@(fes)model_deformation_linear_bbar(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',tet_rule(struct('npts',1)),...
    %         'integration_rule_deviatoric',tet_rule(struct('npts',4))));
    %     eltyd(eix).styl='ch-';
    %     eix=eix+1;
    
    %     eltyd(eix).description ='T10-SRI-bulk';
    %     eltyd(eix).mf =@T10_block;
    %     eltyd(eix).femmf =@(fes)model_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',tet_rule(struct('npts',1)),...
    %         'integration_rule_deviatoric',tet_rule(struct('npts',4)),'split', 'bulk'));
    %     eltyd(eix).styl='kp-';
    %     eix=eix+1;
    %
    %     eltyd(eix).description ='T10-SRI-lambda';
    %     eltyd(eix).mf =@T10_block;
    %     eltyd(eix).femmf =@(fes)model_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',tet_rule(struct('npts',1)),...
    %         'integration_rule_deviatoric',tet_rule(struct('npts',4)),'split', 'lambda'));
    %     eltyd(eix).styl='gp-';
    %     eix=eix+1;
    
    
    
    for eix = 1:length(eltyd)
        
        for mesh= 1:length( meshes )
            
            % Mesh
            [fens,fes]=mesh_import(meshes{mesh},1/1000);
            %         count(fens)
            %                 drawmesh({fens,gcells},'gcells','facecolor','red')
            
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
            nd=abs(mean(ucorn(:,2))/analyt_sol);
            disp([eltyd(eix).description ': ' 'Deflection under the load: ' num2str(nd*100) '%'])
            
            if (graphics )
                clear options
                options.u_scale= u_scale;
                options=deformation_plot_deformation(model_data, options);
            else
                plot(count(fens),nd,eltyd(eix).styl);  hold on
            end
        end
        
    end
end