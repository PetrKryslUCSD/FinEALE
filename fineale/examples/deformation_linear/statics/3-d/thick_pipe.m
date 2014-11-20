% An infinitely long thick walled cylinder with inner boundary radius of 3
% units and outer boundary radius of 9 units is subjected to an internal
% pressure of 1.0 units. A wedge   with unit thickness and a 90-degree
% angle sector is considered for the finite element analysis. The material
% properties are taken as  isotropic linear elastic with E=1000 and
% ?=0·4999 to represent nearly incompressible behavior.
% This problem has been proposed to by MacNeal and Harter [33] as a test on
% an elements ability to
% represent proper behavior for a nearly incompressible material.
% %
% Note: There is an analytical solution to this problem (see the reference below).
%
% Macneal RH, Harder RL (1985) A proposed standard set of problems to test
% finite element accuracy. Finite Elements in Analysis and Design 1: 3-20. 
function thick_pipe
    % Thick-walled cylinder (MacNeal, Harder standard test)
    % Parameters:
    E=1000;
    nu=0.3;
    uzex=  0.004582500000000;
    %     nu=0.49;
    %     uzex= 0.005039925000000;
    %     %     nu=0.499;
    %     %     uzex=  0.005060249250000;
    %         nu=0.4999;
    %         uzex=  0.005062274992500;
    nu=0.499999999;
    uzex= 0.005062499997750;
    R= 6;
    t=6;
    L=t;
    ang=90/180*pi;
    p=  1;
    randshiftmult= 0.37;
    randshiftmult= 0.;
    parshiftmult= 0.5;
    parshiftmult= 0.;
    graphics =  true;
    u_scale=200;
    Plot_stress =  true;
    
    mix = 1;
    clear mesd
    mesd(mix).ref=4;[4:2:10];
    mesd(mix).ny=1;
    mesd(mix).nc=4;
    mesd(mix).nt=1;
    mix = mix+1;
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    %             eltyd(eix).description ='H64';
    %             eltyd(eix).mf =@H64_block;
    %             eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %                 'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
    %             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
    %             eltyd(eix).styl='k*--';
    %             eix=eix+1;
            
                 eltyd(eix).description ='H64-SRI';
            eltyd(eix).mf =@H64_block;
            eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
                'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',3)),...
                'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',4))));
            eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
            eltyd(eix).styl='k*--';
            eix=eix+1;
    %
    %         eltyd(eix).description ='H27';
    %         eltyd(eix).mf =@H27_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
    %         eltyd(eix).styl='ks--';
    %         eix=eix+1;
    
    %     eltyd(eix).description ='T10';% tetrahedron
    %     eltyd(eix).mf =@T10_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
    %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    %     eltyd(eix).styl='k^--';
    %     eix=eix+1;
    
    %     eltyd(eix).description ='T10-SRI';
    %     eltyd(eix).mf =@T10_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',tet_rule(struct('npts',1)),...
    %         'integration_rule_deviatoric',tet_rule(struct('npts',4))));
    %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    %     eltyd(eix).styl='kv-';
    %     eix=eix+1;
    
    %     eltyd(eix).description ='H20R';
    %     eltyd(eix).mf =@H20_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='ro--';
    %     eix=eix+1;
    
    %     %         % Selective reduced integration hexahedron
    %     eltyd(eix).description ='H8-SRI';
    %     eltyd(eix).mf =@H8_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
    %         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='md--';
    %     eix=eix+1;
    
    
    for eix = 1:length(eltyd)
        rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
        
        for mix = 1:length(mesd)
            ns=[];
            uzs=[];
            nfreedofs = [];
            
            for     ref=mesd(mix).ref;
                for     nc=mesd(mix).nc;
                    for     nt=mesd(mix).nt;
                        %% Create the mesh and initialize the geometry
                        [fens,fes]= feval (eltyd(eix).mf, ang, L/2,t,nc,mesd(mix).ny,nt*ref);
                        internal_fenids=fenode_select(fens,struct('box',[0 ang 0 L/2 0 0],'inflate',t/10000));
                        bdry_fes = mesh_boundary(fes, struct('other_dimension',1.0));
                        bcl = fe_select(fens, bdry_fes, ...
                            struct ('box',[0 ang 0 L/2 0 0],'inflate',t/10000));
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
                            % drawmesh({fens,fes}); labels([])
                        end
                        xy=fens.xyz;
                        for i=1:count (fens)
                            a=xy(i,1); y=xy(i,2); r=R-t/2+xy(i,3);
                            xy(i,:)=[r*sin(a) y (r*cos(a))];
                        end
                        fens.xyz=xy;
                        
                        % Compose the model data
                        clear model_data
                        model_data.fens =fens;
                        
                        clear region
                        region.fes= fes;
                        region.femm= eltyd(eix).femmf(fes);
                        model_data.region{1} =region;
                        
                        clear essential
                        essential.component= [1];
                        essential.fixed_value= 0;
                        essential.node_list = fenode_select (fens,struct ('box',[0 0 0 L/2 -100*R 100*R],'inflate',t/10000));
                        model_data.boundary_conditions.essential{1} = essential;
                        
                        clear essential
                        essential.component= [2];
                        essential.fixed_value= 0;
                        essential.node_list = fenode_select (fens,struct ('box',[-100*R 100*R 0 0 -100*R 100*R],'inflate',t/10000));
                        model_data.boundary_conditions.essential{2} = essential;
                        
                        clear essential
                        essential.component= [2];
                        essential.fixed_value= 0;
                        essential.node_list = fenode_select (fens,struct ('box',[-100*R 100*R L/2 L/2 -100*R 100*R],'inflate',t/10000));
                        model_data.boundary_conditions.essential{3} = essential;
                        
                        clear essential
                        essential.component= [3];
                        essential.fixed_value= 0;
                        essential.node_list = fenode_select (fens,struct ('box',[-100*R 100*R 0 L/2 0 0],'inflate',t/10000));
                        model_data.boundary_conditions.essential{4} = essential;
                        
                        clear traction
                        traction.fes =subset(bdry_fes,bcl);;
                        traction.traction= @(x) (p*([1,0,1].*x)'/norm(([1,0,1].*x)));
                        traction.integration_rule =eltyd(eix).surface_integration_rule;
                        model_data.boundary_conditions.traction{1} = traction;
                        
                        % Solve
                        tic; model_data =deformation_linear_statics(model_data);
                        toc
                        
                        %% Transfer the solution to the field u
                        uv=gather_values (model_data.u,internal_fenids);
                        ur=0*internal_fenids;
                        xyz=fens.xyz;
                        for j=1:length(internal_fenids)
                            n=[xyz(internal_fenids(j),1),0,xyz(internal_fenids(j),3)];
                            n=n'/norm(n);
                            ur(j)=uv(j,:)*n;
                        end
                        %                     uzs =[uzs mean(ur)];
                        uzs =[uzs mean(ur)];
                        nfreedofs = [nfreedofs,model_data.u.nfreedofs];
                        
                        %  Exact solution from Felippa''s book
                        a=3; b=9;
                        r=a;
                        ur=p*a^2*(1+nu)*(b^2+r^2*(1-2*nu))/(E*(b^2-a^2)*r)
                        uzex
                        
                        if (graphics )
                            model_data.postprocessing.u_scale= u_scale;
                            model_data.postprocessing.stress_component=3;
                            model_data.postprocessing.stress_range=[-p,+p];
                            model_data=deformation_plot_stress(model_data);
                            %                     options=deformation_plot_deformation(model_data, options);
                        else
                            if (Plot_stress )
                                r  =linspace(a,b,100);
                                idat.component =1; %  Radial component of stress
                                plot(r,p*a.^2/(b^2-a^2).*(1-b^2./r.^2),'k.-','linewidth',3); hold on
                                %                                 idat.component =2; %  Hoop component of stress
                                %                                 plot(r,p*a.^2/(b^2-a^2).*(1+b^2./r.^2),'k.-','linewidth',3); hold on
                                context. output='Cauchy';
                                for ii=1:count(model_data.region{1}.femm.fes)
                                    idat.r =[]; idat.s =[];
                                    idat = inspect_integration_points(model_data.region{1}.femm, ...
                                        model_data.geom, model_data.u, [], ii, context,...
                                        @inspector, idat);
                                    for j  =1:length(idat.r)
                                        plot(idat.r(j),idat.s(j),eltyd(eix).styl,'linewidth',3); hold on
                                    end
                                end
                                title(['Thick pipe, \nu  =' num2str(nu,12)])
                                labels('Radial distance', 'Stress $\sigma_r$')
                                set_graphics_defaults
                            else
                                plot(nfreedofs,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
                            end
                            figure (gcf); pause (1)
                        end
                        
                    end
                end
            end
            disp(['data=[']);
            disp([num2str(nfreedofs)]);
            disp([num2str(uzs)]);
            disp(['];']);
            disp(['h=data(1,:);utip=abs((data(2,:)-uzex)/uzex);']);
            disp(['loglog(h,utip,''' eltyd(eix).styl ''',''lineWIDTH'', 3); hold on;']);
            disp(['Legends{end+1} =''' eltyd(eix).description ''';']);
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