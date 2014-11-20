function axial_wave_expl_1D
    graphics = true;
    igraphics=(1:8:40000000);
    plots = true;
    % Parameters:
    u=physical_units_struct;
    E=205000*u.MEGA*u.PA;
    nu=0.3;
    rho  = 7850*u.KG/u.M^3;
    
    L = 2000*u.MM;
    W = 4*u.MM;
    H = 8*u.MM;
    tolerance  =W/1000;
    vmag = 0.1*u.M/u.SEC;
    tend = 0.0013*u.SEC;
    scale=4000;
    
    %     Mesh
    [fens,fes] = L2_block(L,150,struct('other_dimension',W*H));
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.property = 'isotropic';
    region.E =E;
    region.nu=nu;
    region.rho=rho;
    region.fes= fes;
    region.integration_rule = gauss_rule (struct('dim', 1, 'order', 2));
    model_data.region{1} =region;
    
    clear essential
    essential.component= [];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[L,L],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)0*x+vmag;
    model_data.initial_condition = initial_condition;
    
    model_data.tend = tend;;
    model_data.observer  =@output;;
    corner=fenode_select (fens,struct('box',[0 0],'inflate',tolerance));
    corneru  = [];
    
    % Solve
    model_data =deformation_linear_direct_explicit_CD(model_data);
    
    function output(t, model_data)
        if graphics && (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
            % plot(model_data.geom.values,model_data.u.values);
            % set(gca,'xlim',[0,L]);
            % set(gca,'ylim',[-5e-4,5e-4]/10);
            
            plot(model_data.geom.values,model_data.v.values,'+');
            set_graphics_defaults
            labels('Distance along the bar',' Velocity');
            set(gca,'xlim',[0,L]);
            set(gca,'ylim',[-1,1]*vmag*1.5);
            pause(0.1);
            %             saveas(gcf, [mfilename '-' num2str(etaw) '-' num2str(snapshot) '.png'], 'png');
            igraphics= igraphics(2:end);
        end
    end
end