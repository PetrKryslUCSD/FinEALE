function taylor_bar_expl_nl
    
    
    % Parameters:
    u=physical_units_struct;
    E=117000*u.MEGA*u.PA;
    nu=0.35;
    rho  = 8930*u.KG/u.M^3;
    sigma_y=400*u.MEGA*u.PA; Hi=0*u.MEGA*u.PA;
    iv=227000*u.MM/u.SEC;
    Radius =  3.2*u.MM;
    Length = 32.4*u.MM;
    nperradius=4; nL=28;
    tolerance  =Radius/1000;
    tend = 80e-6*u.SEC;
    scale=4000;
    graphics = true;
    igraphics=10;
    plots = true;
    
    prop = property_deformation_plasticity_linear_hardening(struct('E',E,'nu',nu,'rho',rho,'sigma_y',sigma_y,'Hi',0.0));
    mater = material_deformation_ifr_j2(struct('property',prop));
    
    %     Mesh
    if (1)% Choose hexahedral mesh
       [fens,fes] = H8_quarter_cylinder_n(Radius, Length, nperradius, nL);
        femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
            'integration_rule',gauss_rule(struct('dim',3,'order',2))));
        Surface_integration_rule =gauss_rule(struct('dim',2, 'order',2));
    else% Choose tetrahedral mesh
        [fens,fes] = T10_blockb(L,W,H, 4,1,2);
        [fens,fes] = T10_to_T10MS(fens,fes);
        femm = femm_deformation_nonlinear_t10ms(struct ('material',mater, 'fes',fes, ...
            'integration_rule',tet_rule(struct('npts',1))));
        Surface_integration_rule =tri_rule(struct('npts',1));
    end
    
    %     drawmesh({fens,fes},'fes','facecolor','red'); hold on


    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    
    
    clear region
    region.femm = femm;
    model_data.region{1} =region;
    
    
    clear essential
    essential.component= [1];% symmetry with respect to the YZ plane
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    clear essential
    essential.component= [2];% symmetry with respect to the XZ plane
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[-Inf,Inf,0,0,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{2} = essential;
    clear essential
    essential.component= [3];% clamped end
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[-Inf,Inf,-Inf,Inf,0,0],'inflate',tolerance));
    model_data.boundary_conditions.essential{3} = essential;
    
%% 
% Initial conditions
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)-iv*ones(size(x,1),1)*[0,0,1];
    model_data.initial_condition = initial_condition;
    
    model_data.tend = tend;;
    model_data.observer  =@output;;
    midp=fenode_select (fens,struct('box',[0 0 0 0 0 0],'inflate',tolerance));
    midpu  = [];
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
    nout=0;
    model_data.dt=1e-7;
    tic; model_data =deformation_nonlinear_direct_explicit_CD(model_data);
    toc
    
    function output(t, model_data)
        midpu= [midpu    model_data.un1.reshape(gather_values(model_data.un1, midp))];
        if graphics && (mod(nout,igraphics)==0)
            gv=reset (gv,[]);
            draw(model_data.region{1}.femm.fes, gv, struct ('x', model_data.geom, 'u', model_data.un1,'facecolor','none','edgecolor','blue'));
            view(3);  title(['Time=',num2str(t*1e6)]);
            interact(gv)
            pause(0.1);
        elseif plots && (~graphics)
            figure (pf);
            plot(t/(u.MILLI*u.SEC),midpu(3,end)/u.IN,'rx','Markersize',4); hold on
            labels ('Time [ms]','Deflection [in]')
            pause(0.01)
        end
        nout= nout+1;
        %         end
    end
end