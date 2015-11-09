function sscircularplatepl
disp('Clamped square plate with concentrated force');
% Data listed in the Simo 1990 paper "A class of... "
E=10.92e6;% The Young's modulus is taken  much larger to limit the amount of deflection
nu=0.3;
Magnitude=1.0;
R=5.0;
graphics =  true;
nt=1;
for thickness = 0.1;
    tolerance=0.0001*thickness;
    % Solution listed in the Simo 1990 paper "A class of... "
    analyt_sol= Magnitude/64*((5+nu)/(1+nu) +4/3*(3+nu)/(1-nu^2)*(thickness^2/R^2))*R^4/(E/(1-nu^2)/12*(thickness^3));
    
    
    % Mesh
    Normalized_Deflection =[];
    for nperradius=[2,4,8,12]
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
        
        enl=fenode_select (fens,struct ('box',[0,0, 0,0, 0,thickness],'inflate',tolerance));
        
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
        
        
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        u0=gather_values (model_data.u,enl);
        disp(['Vertical deflection under the load: ' num2str(mean(u0(:,3))) '  --  ' num2str((mean(u0(:,3)))/analyt_sol*100) '%'])
        Normalized_Deflection = [Normalized_Deflection,(mean(u0(:,3)))];
        %         Number_of_elements= [Number_of_elements,n];
        % Plot
        if graphics
            model_data.postprocessing.u_scale= R/8/analyt_sol*(mean(u0(:,3))/analyt_sol);
            %             model_data.postprocessing.show_mesh= 1;
            %             model_data=deformation_plot_deformation(model_data);
                        model_data.postprocessing.output='Cauchy';
                        model_data.postprocessing.stress_component=5;
                        model_data.postprocessing.outputRm=@cylindrical_Rm;
                        model_data.postprocessing.cmap=cadcolors2;
                        %                         model_data=deformation_plot_stress_elementwise(model_data);
                        model_data=deformation_plot_stress(model_data);
        end
    end
    
end

    function Rm= cylindrical_Rm(x,J,label)
        eR=[x(:,1:2),0]'; eZ=[0;0;1];
        eR=eR/norm(eR);
        eT=skewmat(eR)*eZ;
        Rm= [eR,eT,eZ];
    end
end
