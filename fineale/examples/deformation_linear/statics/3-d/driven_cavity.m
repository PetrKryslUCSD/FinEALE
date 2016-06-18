
function driven_cavity
% Driven cavity 
% Parameters:
E=1.0;
   nu =0.4999;
 W= 1; H= 1; T=H/2;
  include_corners = ~true;
graphics = true;
u_scale=0.1;
Plot_stress = true;
Plot_displacement =true;
randshiftmult=0;
      parshiftmult=0;

mix = 1;
clear mesd
mesd(mix).nW=20;
mesd(mix).nH=20;
mesd(mix).nT=1;
mix = mix+1;

prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;

% eltyd(eix).description ='NICE-H8';% tetrahedron
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_nice(struct('fes',fes,'material',mater,...
%     'integration_rule',tensprod_nq_rule(struct('dim',3, 'order',1))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

eltyd(eix).description='H8MSGSO';
eltyd(eix).mf =@H8_block;
eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgso(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2,'order',2));
eix=eix+1;
%
%
% eltyd(eix).description ='H8-Bbar';% tetrahedron
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
%     'pv_bfun',@(p)[1]));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H8-GSRI';% tetrahedron
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H8-GSRI';% tetrahedron
% eltyd(eix).mf =@H8_block_u;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
%     eltyd(eix).description ='H20R';
%     eltyd(eix).mf =@H20_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='ro--';
%     eix=eix+1;

%         eltyd(eix).description ='H64';
%         eltyd(eix).mf =@H64_block;
%         eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%             'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%         eltyd(eix).styl='k*--';
%         eix=eix+1;
%
%         eltyd(eix).description ='H27';
%         eltyd(eix).mf =@H27_block;
%         eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%             'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%         eltyd(eix).styl='ks--';
%         eix=eix+1;

%     eltyd(eix).description ='T10';% tetrahedron
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;

%     eltyd(eix).description ='T10-GSRI';
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%         'integration_rule_constrained',tet_rule(struct('npts',1)),...
%         'integration_rule_unconstrained',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;
%
%     eltyd(eix).description ='H20R';
%     eltyd(eix).mf =@H20_block;
%     eltyd(eix).femmf =@(fes)fem_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='ro--';
%     eix=eix+1;
%
%         eltyd(eix).description ='H27-GSRI-8-27';
%         eltyd(eix).mf =@H27_block;
%         eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%             'integration_rule_constrained', hex_8pt_rule(),...
%             'integration_rule_unconstrained',hex_27pt_rule()...
%             ));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%         eltyd(eix).styl='gd-';
%         eix=eix+1;
%
%         eltyd(eix).description ='H64-GSRI-27-64';
%         eltyd(eix).mf =@H64_block;
%         eltyd(eix).femmf =@(fes)fem_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%             'integration_rule_constrained', hex_27pt_rule(),...
%             'integration_rule_unconstrained',hex_64pt_rule()...
%             ));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%         eltyd(eix).styl='mx-';
%         eix=eix+1;

%     %         % Selective reduced integration hexahedron
%     eltyd(eix).description ='H8-SRI';
%     eltyd(eix).mf =@H8_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
%         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='md--';
%     eix=eix+1;

legends={'Analytical'};
for eix = 1:length(eltyd)
    rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
    
    for mix = 1:length(mesd)
        ns=[];
        uzs=[];
        neqns = [];
        
        for     nW=mesd(mix).nW;
            for     nH=mesd(mix).nH;
                for     nT=mesd(mix).nT;
                    tolerance=min([W,H,T]./[nW,nH,nT])/100;
                    Weps = 2*(1-(include_corners)*1)*tolerance;
                    [fens,fes]= eltyd(eix).mf(W,H,T,nW,nH,nT);
                    if randshiftmult~=0
                        fens = block_random_shift(fens, randshiftmult*[ang/(nc)/3,0,t/(nt*ref)/3],  struct('link', [0, 1, 0]));
                    end
                    if parshiftmult~=0
                        options.displacement=cell(3);
                        ax=parshiftmult*ang/(nc)*2/ref;
                        options.displacement{1} =...
                            @(xyz)(ax*sin(pi*xyz(:,3)/(2*t/(nt*ref))));
                        options.displacement{2} =@(xyz)(0*xyz(:,3));
                        options.displacement{3} =@(xyz)(0*xyz(:,3));
                        fens = block_function_shift(fens, options);
                        % drawmesh({fens,fes},'fes'); labels
                    end
                    
                    % Compose the model data
                    clear model_data
                    model_data.fens =fens;
                    
                    clear region
                    region.fes= fes;
                    region.femm= eltyd(eix).femmf(fes);
                    model_data.region{1} =region;
                    
                    clear essential
                    essential.component= [1,2];
                    essential.fixed_value= 0;
                    essential.node_list = [...
                        fenode_select(fens,struct('box',[0 0 -inf inf -inf inf],'inflate',tolerance)),...
                        fenode_select(fens,struct ('box',[W W -inf inf -inf inf],'inflate',tolerance))];
                    model_data.boundary_conditions.essential{1} = essential;
                    
                    clear essential
                    essential.component= [1,2];
                    essential.fixed_value= 0;
                    essential.node_list = [...
                        fenode_select(fens,struct ('box',[-inf inf 0 0 -inf inf],'inflate',tolerance)),...
                        fenode_select(fens,struct ('box',[-inf inf H H -inf inf],'inflate',tolerance))];
                    model_data.boundary_conditions.essential{2} = essential;
                    
                    clear essential
                    essential.component= [3];
                    essential.fixed_value= 0;
                    essential.node_list = [...
                        fenode_select(fens,struct ('box',[-inf inf -inf inf 0 0],'inflate',tolerance)),...
                    fenode_select(fens,struct ('box',[-inf inf -inf inf T T],'inflate',tolerance))];
                    model_data.boundary_conditions.essential{3} = essential;
                    
                    clear essential
                    essential.component= [1];
                    essential.fixed_value=1.0;
                    essential.node_list = fenode_select(fens,struct('box',[Weps W-Weps H H -inf inf],'inflate',tolerance));
                    model_data.boundary_conditions.essential{4} = essential;
                    
                    % Solve
                    model_data.factorize=false;
                    model_data =deformation_linear_statics(model_data);
                    
                    %                     uzs =[uzs mean(ur)];
                    neqns = [neqns,model_data.u.nfreedofs];
                    
                   
                    
                    if (graphics )
                        model_data.postprocessing.u_scale= u_scale;
                        model_data.postprocessing.output='pressure';
                        model_data.postprocessing.stress_component=1;
                        model_data.postprocessing.stress_range=[-0.4, 0.4];
                        model_data=deformation_plot_stress_elementwise(model_data);
                        %                         model_data=deformation_plot_stress(model_data);
                        % model_data=deformation_plot_deformation(model_data);
                    else
                        if (Plot_stress )
                            fld = field_from_integration_points(model_data.region{1}.femm, model_data.geom, model_data.u, [], [], [], 'pressure', 1);
                            c =fenode_select (fens,struct('box',[-inf,inf,H/2, H/2,-inf,inf],'inflate',tolerance));
                            pressurex=fld.values(c);
                            positionx=model_data.geom.values(c,1);
                            [ignore,ix]=sort(positionx(:,1));
                            px=positionx(ix,1);
                            press=pressurex(ix);
                            plot(px,press,name_to_style(eltyd(eix).description),'linewidth',3); hold on
                            labels(  'Distance', 'p')
                        elseif (Plot_displacement)
                            c =fenode_select (fens,struct('box',[W/2, W/2,-inf,inf,-inf,inf],'inflate',tolerance));
                            ux=model_data.u.values(c,1);
                            positionx=model_data.geom.values(c,2);
                            [ignore,ix]=sort(positionx(:,1));
                            px=positionx(ix,1);
                            displacement=ux(ix);
                            plot(px,displacement,name_to_style(eltyd(eix).description),'linewidth',3); hold on
                            c =fenode_select (fens,struct('box',[-inf,inf,H/2,H/2,-inf,inf],'inflate',tolerance));
                            uy=model_data.u.values(c,2);
                            positionx=model_data.geom.values(c,1);
                            [ignore,ix]=sort(positionx(:,1));
                            px=positionx(ix,1);
                            displacement=uy(ix);
                            plot(px,displacement,name_to_style(eltyd(eix).description),'linewidth',3); hold on
                            grid on
                            labels(  'Distance', 'u_x')
                        else
                            plot(neqns,uzs/uzex,name_to_style(eltyd(eix).description),'linewidth',3); hold on
                        end
                        figure (gcf); pause (1)
                    end
                    
                end
            end
        end
        disp(['data=[']);
        disp([num2str(neqns)]);
        disp([num2str(uzs)]);
        disp(['];']);
        disp(['h=data(1,:);utip=abs((data(2,:)-uzex)/uzex);']);
        disp(['loglog(h,utip,''' name_to_style(eltyd(eix).description) ''',''lineWIDTH'', 3); hold on;']);
        disp(['Legends{end+1} =''' eltyd(eix).description ''';']);
        legends{end+1} =eltyd(eix).description;
    end
    if (~graphics)
        %         legend (legends);
        
    end
end

    function idat =inspector(idat, out, xyz, pc)
        %         if (isempty(idat.r))
        Projection=[xyz(1),0,xyz(3)];
        r=norm(Projection);
        e1p=Projection'/r;
        e3p=[0,1,0]';
        e2p=skewmat(e3p)*e1p;
        Rm= [e1p,e2p,e3p];
        tm = stress_6v_to_3x3t (mater,out);
        tpm = Rm'*tm*Rm;
        sp = stress_3x3t_to_6v (mater,tpm);
        idat.r(end+1) =r;
        idat.s(end+1)=sp(idat.component);
        %         end
    end
end