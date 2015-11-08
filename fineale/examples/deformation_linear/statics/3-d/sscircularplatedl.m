function sscircularplatedl
disp('Simply supported square plate with distributed load');
% Data listed in the Simo 1990 paper "A class of... "
    E=10.92;
nu=0.3;
Magnitude=1;
R=5;
graphics =  false;
for thickness = 0.1;
    tolerance=0.0001*thickness;
    % analytical solution for the vertical deflection under the load
    % Zienkiewicz, Taylor, The finite element method, fifth edition, volume 2,
    %
    %     analyt_sol= 1/64*((5+nu)/(1+nu))/(E/(1-nu^2)/12*(thickness^3));
    % Solution listed in the Simo 1990 paper "A class of... "
    analyt_sol= Magnitude/64*((5+nu)/(1+nu) +4/3*(3+nu)/(1-nu^2)*(thickness^2/R^2))*R^4/(E/(1-nu^2)/12*(thickness^3));
    
    nt=1;
    % Mesh
    Normalized_Deflection =[];
    for nperradius=[4,8,16,32]
        nt=nt+1;
        [fens,fes] = Q4_circle_n(R, nperradius, 1.0);
        [fens,fes] = H8_extrude_Q4(fens,fes,nt,@(x,k)([x(1),x(2),k*thickness/nt]));
        count(fes)/nt
        prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
        mater = material_deformation_linear_triax (struct('property',prop ));
        femm=femm_deformation_nonlinear_h8msgso(struct('fes',fes,'material',mater,...
            'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
        bdry_fes = mesh_boundary(fes);
        topl=fe_select(fens, bdry_fes, ...
            struct ('box',[-inf inf -inf inf  thickness thickness],'inflate',tolerance));
        botl=fe_select(fens, bdry_fes, ...
            struct ('box',[-inf inf -inf inf  0 0],'inflate',tolerance));
        x0l=fe_select(fens, bdry_fes, ...
            struct ('box',[0 0 -inf inf  0 thickness],'inflate',tolerance));
        y0l=fe_select(fens, bdry_fes, ...
            struct ('box',[-inf inf 0 0   0 thickness],'inflate',tolerance));
        cyll=setdiff(1:count(bdry_fes),[topl,botl,x0l,y0l]);
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.fes= fes;
        region.femm= femm;
        model_data.region{1} =region;
        
        clear essential
        essential.component= [1];
        essential.fixed_value= 0;
        essential.node_list = connected_nodes (subset(bdry_fes, x0l));
        model_data.boundary_conditions.essential{1} = essential;
        
        clear essential
        essential.component= [2];
        essential.fixed_value= 0;
        essential.node_list = connected_nodes (subset(bdry_fes, y0l));
        model_data.boundary_conditions.essential{2} = essential;
        
        clear essential
        essential.component= [3];
        essential.fixed_value= 0;
        essential.node_list = connected_nodes (subset(bdry_fes, cyll));
        model_data.boundary_conditions.essential{3} = essential;
        
        
        clear traction
        traction.fes =subset(bdry_fes,topl);;
        traction.traction= [0;0;Magnitude];
        traction.integration_rule =gauss_rule(struct('dim',2, 'order',2));
        model_data.boundary_conditions.traction{1} = traction;
        
        
        enl=fenode_select (fens,struct ('box',[0,0, 0,0, 0,thickness],'inflate',tolerance));
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        u0=gather_values (model_data.u,enl);
        disp(['Vertical deflection under the load: ' num2str(mean(u0(:,3))) '  --  ' num2str((mean(u0(:,3)))/analyt_sol*100) '%'])
        Normalized_Deflection = [Normalized_Deflection,(mean(u0(:,3)))];
        %         Number_of_elements= [Number_of_elements,n];
        % Plot
        if graphics
            model_data.postprocessing.u_scale= R/4/(mean(u0(:,3)));
                model_data.postprocessing.show_mesh= 1;
                model_data=deformation_plot_deformation(model_data);
        end
    end
    Normalized_Deflection
    Normalized_Deflection/analyt_sol*100
end