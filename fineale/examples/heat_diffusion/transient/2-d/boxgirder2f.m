function boxgirder2f
% Transient temperature field on a box go to bridge during a five day period.
kappa_concrete=1.81*[1 0; 0 1]; % W/K/m
rho_concrete = 2350;% kg/m^3
cv_concrete =0.22*rho_concrete;
Ti = 19;% in degree Celsius
Ta_inner=[19, 18, 17, 18, 20, 26, 30, 30, 29, 26, 22, 20];
Ta_inner= [Ta_inner Ta_inner Ta_inner Ta_inner Ta_inner];
Ta_lower=[19, 18, 17, 25, 28, 35, 32, 31, 32, 27, 22, 20];
Ta_lower= [Ta_lower Ta_lower Ta_lower Ta_lower Ta_lower];
Ta_upper=[19, 18, 17, 35, 65, 80, 80, 60, 45, 35, 22, 20];
Ta_upper= [Ta_upper Ta_upper Ta_upper Ta_upper Ta_upper];
Ta_vert=[19, 18, 17, 52, 65, 56, 35, 31, 32, 27, 22, 20];
Ta_vert= [Ta_vert Ta_vert Ta_vert Ta_vert Ta_vert];
hs=[6.2, 6.2, 6, 6.3, 7.0, 7.3, 9.5, 13.4, 17.1, 17.6, 15.4, 9.5];
hs= [hs hs hs hs hs];
%
times =(0:2:2*(length(hs)-1))*3600 ; % in seconds;
dt=2*3600; % time  increment in seconds
tend= (4*24+16)*3600; % length of the time interval
t=0;
crank.hot = 'co-'; crank.cold ='m*-'; crank.theta = 0.5;
bweul.hot = 'bo-'; bweul.cold ='r*-'; bweul.theta = 1.0;
settings =bweul;
theta = settings.theta; 
online_graphics=true;% plot the solution as it is computed?
mesh_size=0.1; %
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    'curve 1 line 0 0 2.60 0',...
    'curve 2 line 2.60 0 2.6 2.00',...
    'curve 3 line 2.6 2.0 4.85 2.0',...
    'curve 4 line 4.85 2 4.85 2.2',...
    'curve 5 line 4.85 2.2 0.0 2.2',...
    'curve 6 line 0.0 2.2 -4.85 2.2 ',...
    'curve 7 line -4.85 2.2 -4.85 2',...
    'curve 8 line -4.85 2.0 -2.6 2.0 ',...
    'curve 9 line -2.6 2.00 -2.60 0 ',...
    'curve 10 line -2.60 0 0 0',...
    'curve 11 line -2.1 2 2.1 2',...
    'curve 12 line 2.1 2 2.1 0.2',...
    'curve 13 line 2.1 0.2 -2.1 0.2',...
    'curve 14 line -2.1 0.2 -2.1 2',...
    ['subregion 1  property 1 boundary '...
    ' 1 2 3 4 5 6 7 8 9 10 hole 11 12 13 14'],...
    ['m-ctl-point constant ' num2str(mesh_size)]
    }, 1.0);
% drawmesh({fens,fes},'fes')
prop_concrete=...
    property_heat_diffusion (struct('thermal_conductivity',kappa_concrete,...
    'specific_heat',cv_concrete,'source',0.0));
mater_concrete=material_heat_diffusion (struct('property',prop_concrete));
feb_concrete = femm_heat_diffusion (struct ('material',mater_concrete,...
    'fes',subset(fes,groups{1}),...
    'integration_rule',tri_rule(struct('npts',1))));
efes_inner=subset(edge_fes,[edge_groups{[(11:14)]}]);
efes_lower=subset(edge_fes,[edge_groups{[1, 3, 7:10]}]);
efes_upper=subset(edge_fes,[edge_groups{[5, 6]}]);
efes_vert=subset(edge_fes,[edge_groups{[2, 4]}]);
geom = nodal_field(struct('name',['geom'], 'dim', 2, 'fens',fens));
tempn=nodal_field(struct('name',['tempn'], 'dim', 1, 'nfens',...
    geom.nfens));
amb = clone(tempn, ['amb']);

tempn = numberdofs (tempn);
tempn = scatter_sysvec(tempn,gather_sysvec(tempn)*0+Ti);
Cm = capacity(feb_concrete, sysmat_assembler_sparse, geom, tempn);
Kmc = conductivity(feb_concrete, sysmat_assembler_sparse, geom, tempn);

if online_graphics
    gv=graphic_viewer;
    gv=reset (gv,struct('limits',[-4.85, 4.85, 0, 2.2, 0, 2.5]));%,'peek',true));
    dcm=data_colormap(struct('range',[20, 70],'colormap',jet));
    frame=1;
end
minmaxT = []; ts = [];
t= 0; 
while t<tend+0.1*dt % Time stepping
    if online_graphics
        gv=reset (gv,struct('limits',[-4.85, 4.85, 0, 2.2, 0, 3.5]));%,'peek',~true));
        set(gca,'FontSize', 14)
        T=tempn.values;
        colorfield=nodal_field(struct('name',['colorfield'],'data',...
            map_data(dcm, T)));
        geomT=nodal_field(struct ('name', ['geomT'], ...
            'data',[geom.values, 0.1*(tempn.values-Ti)]));
        draw(fes,gv,struct('x',geomT, 'u',0*geomT,...
                'colorfield',colorfield, 'edgecolor','none'));
            draw(fes,gv,struct('x',geom, 'u',0*geom, ...
                'facecolor','none'));
        camset(gv,[27.6991   19.1811   28.9747    0.1788    1.1224    2.8844   -0.3120    0.9029   -0.2958 9.0864])
        xlabel('X [m]');        ylabel('Y [m]');
        zl=zlabel('Temperature [{}^0{C}]');
        set(zl,'Rotation',0);
        hour =mod (round(t/3600), 24);
        day =floor(t/3600/24)+1;
        title (['Day ' num2str(day) ', ' num2str(hour) ' hours']); hold off;
        set(gcf,'renderer','zbuffer');pause(1)
        gif_animation_add_frame(gcf,frame,'Moviefilename.gif',340) 
        frame = frame+1;
        
        %         saveas(gcf, ['shrinkfit-' num2str(t) '.png'], 'png');
    end
    tempn1 = tempn;
    h=interp1(times,hs,t);
    efeb = femm_heat_diffusion (struct ('material',mater_concrete,...
        'fes',subset(edge_fes,[edge_groups{:}]),...
        'integration_rule',gauss_rule(struct('dim',1,'order',1)),...
        'surface_transfer', h));
    Km = Kmc+surface_transfer(efeb, sysmat_assembler_sparse, geom, tempn);
    conn = connected_nodes(efes_inner);
    Ta = interp1(times,Ta_inner,t);
    amb = set_ebc(amb, conn, conn*0+1, [], conn*0+Ta);
    conn = connected_nodes(efes_lower);
    Ta = interp1(times,Ta_lower,t);
    amb = set_ebc(amb, conn, conn*0+1, [], conn*0+Ta);
    conn = connected_nodes(efes_upper);
    Ta = interp1(times,Ta_upper,t);
    amb = set_ebc(amb, conn, conn*0+1, [], conn*0+Ta);
    conn = connected_nodes(efes_vert);
    Ta = interp1(times,Ta_vert,t);
    amb = set_ebc(amb, conn, conn*0+1, [], conn*0+Ta);
    amb = apply_ebc (amb);
    F= surface_transfer_loads(efeb,sysvec_assembler, geom,tempn,amb);
    Tn=gather_sysvec(tempn);
    minmaxT= [minmaxT; min(Tn),max(Tn)];
    ts= [ts,t];
    Tn1=(1/dt*Cm+theta*Km)\((1/dt*Cm-(1-theta)*Km)*Tn+F);
    tempn = scatter_sysvec(tempn1,Tn1);
    t=t+dt;
%     if t>24*3*tmul
%         dt=0.1*tmul; % time 0
%     end
end
if online_graphics
    interact(gv)
    figure(gcf)
end
plot(ts/3600,minmaxT(:, 1),settings.hot,'linewidth', 2); hold on
plot(ts/3600,minmaxT(:, 2),settings.cold,'linewidth', 2)
set(gca,'FontSize', 14); grid on
xlabel('Time in hours')
ylabel('Temperature in degrees Celsius')
