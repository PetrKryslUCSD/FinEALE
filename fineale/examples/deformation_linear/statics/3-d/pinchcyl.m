% The pinched  cylinder benchmark

% Macneal RH, Harder RL (1985) A proposed standard set of problems to test
% finite element accuracy. Finite Elements in Analysis and Design 1: 3-20. 
function gv =pinchcyl
    E=3e6;
    nu=0.3;
    thickness = 3.0;
    % analytical solution for the vertical deflection under the load
    uzex=-1.82488e-5;
    graphics = true;
    % Mesh
    R=300;
    L= 600;
    nt=1;
    tolerance =thickness/nt/10;
    load = [0,0,1]';
    randshiftmult=0;
    parshiftmult=0.;
    u_scale=10000000;
    
    mix = 1;
    clear mesd
    mesd(mix).ref=[2:2:10];
    mesd(mix).nl=1;
    mesd(mix).nt=1;
    mix = mix+1;
    
    clc;
    disp(['% ' mfilename]);
    disp(['nt=' num2str(nt,12) ';']);
    disp(['parshiftmult=' num2str(parshiftmult,12) ';']);
    disp(['uzex=' num2str(uzex,12) ';']);
    
    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    eltyd(eix).description ='H20R';
    eltyd(eix).mf =@H20_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).volume_integration_rule=gauss_rule(struct('dim',3, 'order',2));
    eltyd(eix).styl='ro--';
    eix=eix+1;
    %
    %     eltyd(eix).description ='H20';
    %     eltyd(eix).mf =@H20_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='r.--';
    %     eix=eix+1;
    %
    %     eltyd(eix).description ='H27';
    %     eltyd(eix).mf =@H27_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='rh--';
    %     eix=eix+1;
    
    
    %
    mix=1;
    
    legends ={};
    for eix = 1:length(eltyd)
        rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
        ns=[];
        nfens=[];
        uzs =[];
        for n  = mesd(mix).ref
            nl=mesd(mix).nl;
            nt=mesd(mix).nt;
            % Mesh
            [fens,fes]= eltyd(eix).mf(90/360*2*pi,L/2,thickness,n,n,nt);
            % Shape into a cylinder
            if randshiftmult~=0
                fens = block_random_shift(fens, randshiftmult*[ang/(nl*n)/3,R/(nl*n)/3, 0],  struct('link', [0, 0, 1]));
            end
            if parshiftmult~=0
                ang =90/360*2*pi;
                options.displacement=cell(3);
                ax=0.5*parshiftmult*ang/(nl*n)/2;
                options.displacement{1} =@(xyz)(ax*cos(1.2*pi*xyz(:,2)/(L/(nl*n))));
                ay=0.4*parshiftmult*L/(nl*n);
                options.displacement{2} =@(xyz)(ay*sin(0.6*pi*xyz(:,1)/(ang/(nl*n))));
                options.displacement{3} =@(xyz)(0*xyz(:,3));
                fens = block_function_shift(fens, options);
            end
            xyz=fens.xyz;
            for i=1:count (fens)
                a=xyz(i,1); y=xyz(i,2); z=xyz(i,3);
                xyz(i,:)=[(R-thickness/2+z)*sin(a) y (R-thickness/2+z)*cos(a)];
            end
            fens.xyz=xyz;
            %         count(fens)
            %                                        gv= drawmesh({fens,fes},'fes','facecolor','red');;
            %             rotate(gv); figure(gcf)
            % Compose the model data
            clear model_data
            model_data.fens =fens;
            
            clear region
            region.fes= fes;
            region.femm= eltyd(eix).femmf(fes);
            model_data.region{1} =region;
            
            clear essential
            essential.component= [1,3];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[-10000 10000 0 0 -10000 10000],'inflate',tolerance));
            model_data.boundary_conditions.essential{1} = essential;
            
            clear essential
            essential.component= [2];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[-10000 10000 L/2 L/2 -10000 10000],'inflate',tolerance));
            model_data.boundary_conditions.essential{2} = essential;
            
            clear essential
            essential.component= [1];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[0 0 -10000 10000 -10000 10000],'inflate',tolerance));
            model_data.boundary_conditions.essential{3} = essential;
            
            clear essential
            essential.component= [3];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[-10000 10000 -10000 10000 0 0],'inflate',tolerance));
            model_data.boundary_conditions.essential{4} = essential;
             
            clear nodal_force
            nodal_force.node_list =fenode_select (fens,struct ('box',[0 0 L/2 L/2 -1000  1000],'inflate',tolerance));
            nodal_force.force= [0;0;-1/4/length(nodal_force.node_list)];
            model_data.boundary_conditions.nodal_force{1} = nodal_force;
            
            % Solve
            model_data =deformation_linear_statics(model_data);
            
            corn=fenode_select (fens,struct('box',[0 0 L/2 L/2 -10000 10000],'inflate',tolerance));
            uv= gather_values(model_data.u, corn);
            uz=sum(uv(:,3))/length(corn) ;
            nd=abs(uz/uzex);
            disp(['% ' eltyd(eix).description ': ' 'Deflection under the load: ' num2str(nd*100) '%'])
            uzs =[uzs , uz]; ns =[ns , n]; 
                        nfens = [nfens,model_data.u.nfens];
            
            if (graphics )
                model_data.postprocessing.u_scale= u_scale;
                model_data=deformation_plot_deformation(model_data);
            else
                
            end
        end
        if (graphics )
        else
            %             semilogx(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
            plot(ns,abs(uzs/uzex),eltyd(eix).styl,'linewidth',3);  hold on
            title(['Pinched cylindrical shell, nt  =' num2str(nt)])
            figure (gcf); pause (1)
        end
        legends{end+1} =eltyd(eix).description;
        format long
        disp(['Data{end+1}=[']);
        disp([num2str(ns,15)]);
        disp([num2str(uzs,15)]);
        disp([num2str(nfens,15)]);
        disp(['];' 'Style{end+1}=''' eltyd(eix).styl ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
    end
    if (~graphics )
        legend (legends);
        labels(  'Number of equations', 'Estimated true error')
        grid on
        set_graphics_defaults
    end
end