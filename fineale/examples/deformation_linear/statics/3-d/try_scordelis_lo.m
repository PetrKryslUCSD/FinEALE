% As a test case a cylindrical shell roof was chosen as here both membrane
% and bending stresses  are of importance. The particular example
% originally solved by Scordelis and co-workers  has been
% used frequently since for assessment of finite element performance
% and hence will be useful here for comparison purposes.
%
% Reference:
% Zienkiewicz OC, Taylor RL, Too JM (1971) Reduced integration technique in
% general analysis of plates and shells. International Journal for
% Numerical Methods in Engineering 3: 275-290.
% Original reference:
% A. C. Scordelis and K. S. Lo, ‘Computer analysis of cylindrical shells’, ACI Jnl, 61, 539-561 (1964).
%
function gv =scordelis_lo
    E=4.32e8;
    nu=0;
    R=25;
    L=2*R;
    t= 0.25;
    tolerance =t/1e6;
    % t= 0.5;
    ang=40/180*pi;
    p=  90/t;
    load=[0;0;-p];
    uzex=-0.3024;%  This is in feet, in inches: 3.6288
    randshiftmult= 0.27;
    randshiftmult= 0.;
    % parshiftmult=0.2;
    parshiftmult=0;
    graphics = true;
    u_scale=20;
    
    mix = 1;
    clear mesd
    mesd(mix).ref=[3];
    mesd(mix).nl=1;
    mesd(mix).nt=1;
    mix = mix+1;
    
    clc;
    disp(['% ' mfilename]);
    disp(['uzex=' num2str(uzex,12) ';']);
    
    % Material
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    %         eltyd(eix).description ='H20R';
    %         eltyd(eix).mf =@H20_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %         eltyd(eix).volume_integration_rule=gauss_rule(struct('dim',3, 'order',2));
    %         eix=eix+1;
    %
    %         eltyd(eix).description ='H20';
    %         eltyd(eix).mf =@H20_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %         eltyd(eix).volume_integration_rule=gauss_rule(struct('dim',3, 'order',3));
    %         eltyd(eix).styl='r.--';
    %         eix=eix+1;
    
    %     eltyd(eix).description ='H27';
    %     eltyd(eix).mf =@H27_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='rh--';
    %     eix=eix+1;
    
                eltyd(eix).description ='T10';
                eltyd(eix).mf =@T10_blockca;
                eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
                    'integration_rule',tet_rule(struct('npts',4))));
                eltyd(eix).volume_integration_rule=tet_rule(struct('npts',4));
                eix=eix+1;
    
    %             eltyd(eix).description ='T4';
    %             eltyd(eix).mf =@T4_blockb;
    %             eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %                 'integration_rule',tet_rule(struct('npts',1))));
    %             eltyd(eix).volume_integration_rule=tet_rule(struct('npts',1));
    %             eix=eix+1;
    
    %
    mix=1;
    
    legends ={};
    for eix = 1:length(eltyd)
        rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
        ns=[];
        uzs =[];
        for n  = mesd(mix).ref
            nl=mesd(mix).nl;
            nt=mesd(mix).nt;
            % Mesh
            [fens,fes]= eltyd(eix).mf(ang, L/2,t,nl*n,nl*n,nt*n);
            free_edge_fenids=fenode_select(fens,struct('box',[ang, ang, L/2 L/2 -100*R 100*R],'inflate',tolerance));
            if randshiftmult~=0
                fens = block_random_shift(fens, randshiftmult*[ang/(nl*n)/3,R/(nl*n)/3, 0],  struct('link', [0, 0, 1]));
            end
            if parshiftmult~=0
                options.displacement=cell(3);
                ax=parshiftmult*ang/(nl*n);
                options.displacement{1} =@(xyz)(ax*sin(2*pi*xyz(:,2)/(ang/(nl*n))));
                ay=parshiftmult*L/2/(nl*n);
                options.displacement{2} =@(xyz)(ay*sin(19*pi*xyz(:,1)/(L/2/(nl*n))));
                options.displacement{3} =@(xyz)(0*xyz(:,3));
                fens = block_function_shift(fens, options);
                %                         mesh{1}=fens;
                %                         mesh{2}=fes;
                %                         drawmesh(mesh); view(3); pause%(1)
            end
            xy=fens.xyz;
            for i=1:count (fens)
                a=xy(i,1); y=xy(i,2); r=R-t/2+xy(i,3);
                xy(i,:)=[r*sin(a) y (r*cos(a)-R)];
            end
            fens.xyz=xy;
            %         count(fens)
            %                             drawmesh({fens,fes},'fes','facecolor','red')
            
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
            essential.node_list = fenode_select (fens,struct ('box',[0 100*R 0 0 -100*R 100*R],'inflate',tolerance));
            model_data.boundary_conditions.essential{1} = essential;
            
            clear essential
            essential.component= [2];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[0 100*R L/2 L/2 -100*R 100*R],'inflate',tolerance));
            model_data.boundary_conditions.essential{2} = essential;
            
            clear essential
            essential.component= [1];
            essential.fixed_value= 0;
            essential.node_list = fenode_select (fens,struct ('box',[0 0 0 L/2 -100*R 100*R],'inflate',tolerance));
            model_data.boundary_conditions.essential{3} = essential;
            
            clear body_load
            body_load.force=load;
            body_load.fes  =fes;;
            body_load.integration_rule  = eltyd(eix).volume_integration_rule;
            model_data.body_load{1} = body_load;
            
            % Solve
            model_data =deformation_linear_statics(model_data);
            
            uv= gather_values(model_data.u, free_edge_fenids);
            uz=sum(uv(:,3))/length(free_edge_fenids) ;
            nd=abs(uz/uzex);
            disp(['% ' eltyd(eix).description ': ' 'Deflection under the load: ' num2str(nd*100) '%'])
            uzs =[uzs , uz]; ns =[ns , n];
%             nfreedofs = [nfreedofs,model_data.u.nfreedofs];
                    
            if (graphics )
                model_data.postprocessing.u_scale= u_scale;
                model_data.postprocessing.camera  = 1.0e+02 *[ 0.815002706544087   1.472765509480047   0.374529960567960   0.056880983027792   0.131188079971008 -0.027124953925793  -0.001242969037269  -0.002194705642988   0.009676688230641   0.103395849072020];
                model_data=deformation_plot_deformation(model_data);
            else
                
            end
        end
        if (graphics )
        else
            %             semilogx(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
            plot(ns,abs(uzs/uzex),name_to_style(eltyd(eix).description),'linewidth',3);  hold on
            %             title(['Scordelis-Lo'])
            figure (gcf); pause (1)
        end
        legends{end+1} =eltyd(eix).description;
        format long
        disp(['Data{end+1}=[']);
        disp([num2str(ns,15)]);
        disp([num2str(uzs,15)]);
        disp(['];' 'Style{end+1}=''' name_to_style(eltyd(eix).description) ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
    end
    if (~graphics )
        legend (legends);
    labels(  'Number of elements per edge', 'Norm. deflection')
    grid on
    end
    set_pub_defaults
end