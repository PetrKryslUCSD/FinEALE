function bending_wave_expl_nl
    % Bending wave in a clamped beam. Centered-difference integration.
    graphics = false;
    igraphics=(1:100:40000000);
    plots = true;
    % Parameters:
    u=physical_units_struct;
    E=205000*u.MEGA*u.PA;
    nu=0.3;
    rho  = 7850*u.KG/u.M^3;
    eta=2000*u.PA*u.SEC;
    loss_tangent=0.0001;
    G=E/2/(1+nu);
    frequency = 1/0.0058;
    eta = loss_tangent*G/(2*pi*frequency);
    L = 200*u.MM;
    W = 4*u.MM;
    H = 8*u.MM;
    tolerance  =W/1000;
    vmag = 0.1*u.M/u.SEC;
    tend = 0.0013*u.SEC;
    scale=4000;
    
    %     Mesh
    [fens,fes] = H8_block(L,W,H, 50,1,4);
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    prop = property_deformation_linear_iso (struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_stvk_triax(struct('property',prop));
    %     femm = femm_deformation_nonlinear_h8msgs(struct ('material',mater, 'fes',fes, ...
    %         'integration_rule',gauss_rule(struct('dim',3,'order',2))));
femm = femm_deformation_nonlinear(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));

    clear region
    region.femm = femm;
    model_data.region{1} =region;
    
    
    clear essential
    essential.component= [];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[L,L,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)0*x+ones(size(x,1),1)*[0,0,vmag];
    model_data.initial_condition = initial_condition;
    
    model_data.tend = tend;;
    model_data.observer  =@output;;
    corner=fenode_select (fens,struct('box',[0 0 0 0 0 0],'inflate',tolerance));
    corneru  = [];
    if graphics
        gv=graphic_viewer;
        gv=reset (gv,[]);
        peekb  = uicontrol('Parent',gcf,'Style','pushbutton',...
                    'Units','Normalized','Position',[.9 .0 .1 .05],...
                    'String','Peek','Value',0,...
                    'TooltipString','Invoke command line in order to peek at data',...
                    'Callback',@(h,varargin)keyboard);
    elseif plots
        pf = figure(gcf);
    end
    
    % Solve
%         model_data.dt=1e-8;
    tic; model_data =deformation_nonlinear_direct_explicit_CD(model_data);
    toc
    
    function output(t, model_data)
        disp(num2str(t))
        if graphics && (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
            %         clear(gv);
            gv=reset (gv,[]);
            draw(model_data.region{1}.femm.fes, gv, struct ('x', model_data.geom, 'u', scale*model_data.u,'facecolor','none','edgecolor','blue'));
            view(3); pause(0.1);
            %             saveas(gcf, [mfilename '-' num2str(etaw) '-' num2str(snapshot) '.png'], 'png');
            igraphics= igraphics(2:end);
        elseif plots
            corneru= [corneru    model_data.u.reshape(gather_values(model_data.u, corner))];
            if (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
                figure (pf);
                %             plot (corneru(1,:),corneru(2,:),'g'); hold on
                plot(t,corneru(3,end),'ro','Markersize',4); hold on
                igraphics= igraphics(2:end);
            end
        end
    end
end