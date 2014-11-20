% The curved cantilever beam with two different end loading conditions is
%  studied. The beam is loaded in-plane or out-of-plane.
%
% Reference:
%   R.  H. MacNeal  and  R.  L.  Harder, ‘A  proposed  standard set  of  problems  to test  finite element accuracy’,  Finite
% Elements  Anal. Des., 1, 3-20  (1985).

function curved_beam
    E=1.0e7;
    nu=0.25;
    rin=4.12;
    rex =4.32;
    t= 0.1;
    p=  1.0/(rex-rin)/t;
    load =[0;0;0];% out of plane
    udex = 0.498301636614788;% Harder: .5022;% out of plane
    dir=3;% out of plane
    %         udex = 0.088519352226485;% Harder:.08734;% in plane
    %         dir=2;% in plane
    load(dir)=p;
    htol=min([t])/1000;
    graphics = true;
    stress_range = [-30,30];
    u_scale=10;
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    format long
    clc
    disp(['% ' mfilename]);
    disp(['dir=' num2str(dir,12) ';']);
    disp(['udex=' num2str(udex,12) ';']);
    
    clear mesd
    mix = 1;
    mesd(mix).ref= [1,2,3];%[1,2,4];
    mesd(mix).nl=6;
    mesd(mix).nt=1;
    mix = mix+1;
    
    clear eltyd
    eix=1;
    
    eltyd(eix).description ='H20R';
    eltyd(eix).mf =@H20_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
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
    
    legends ={};
    
    mix = 1;
    for eix = 1:length(eltyd)
        
        ns=[];
        uzs=[];
        nfreedofs = [];
        
        for     ref=mesd(mix).ref;
            for     nl=mesd(mix).nl;
                for     nt=mesd(mix).nt;
                    %% Create the mesh and initialize the geometry
                    [fens,fes]= eltyd(eix).mf(rex-rin,pi/2,t,nt*ref,nl*ref,nt*ref);
                    bdry_fes = mesh_boundary(fes, []);
                    icl = fe_select(fens, bdry_fes, struct('box', [-100*rex, 100*rex,pi/2,pi/2,-100*rex, 100*rex],'inflate',htol));
                    xy=fens.xyz;
                    for i=1: size(xy,1)
                        r=rin+xy(i,1); a=xy(i,2);
                        xy(i,:)=[r*cos(a) r*sin(a) xy(i,3)];
                    end
                    fens.xyz=xy;
                    %                         drawmesh({fens,fes},'shrink', 0.8,'facecolor','red');
                    
                    % Compose the model data
                    clear model_data
                    model_data.fens =fens;
                    
                    clear region
                    region.fes= fes;
                    region.femm= eltyd(eix).femmf(fes);
                    model_data.region{1} =region;
                    
                    clear essential
                    essential.component= [1,2,3];
                    essential.fixed_value= 0;
                    essential.node_list = fenode_select(fens,struct('box',[-100*rex 100*rex 0 0 -100*rex 100*rex],'inflate',htol));
                    model_data.boundary_conditions.essential{1} = essential;
                    
                    clear traction
                    bdry_fes = mesh_boundary(fes, []);
                    traction.fes =subset(bdry_fes,icl);;
                    traction.traction= load;
                    traction.integration_rule =eltyd(eix).surface_integration_rule;
                    model_data.boundary_conditions.traction{1} = traction;
                    
                    % Solve
                    model_data =deformation_linear_statics(model_data);
                    
                    
                    free_end_fenids=fenode_select(fens,struct('box',[0 0 -100*rex 100*rex -100*rex 100*rex],'inflate',htol));
                    uv=gather_values (model_data.u,free_end_fenids);
                    uz=mean (abs(uv(:,dir)));
                    disp (['% ' eltyd(eix).description  ': u =' num2str(uz) ', ' num2str(uz/udex*100) ' %'])
                    uzs =[uzs uz];
                    nfreedofs = [nfreedofs,model_data.u.nfreedofs];
                    
                    if (graphics )
                        model_data.postprocessing.u_scale= u_scale;
                        model_data.postprocessing.stress_component=2;
                        %                     options=deformation_plot_stress(model_data, options)
                        model_data=deformation_plot_deformation(model_data);
                    else
                        %                     plot(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
                        %                     figure (gcf); pause (1)
                    end
                end
            end
        end
        if (graphics )
        else
            %             semilogx(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
            loglog(nfreedofs,abs((uzs-udex)/udex),eltyd(eix).styl,'linewidth',3); hold on
            title(['curved beam, dir=' num2str(dir)])
            figure (gcf); pause (1)
        end
        legends{end+1} =eltyd(eix).description;
        format long
        disp(['Data{end+1}=[']);
        disp([num2str(nfreedofs,15)]);
        disp([num2str(uzs,15)]);
        disp(['];' 'Style{end+1}=''' eltyd(eix).styl ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
        legend (legends);
        labels(  'Number of equations', 'Estimated true error')
        grid on
        set_graphics_defaults
    end
end