function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.932;
convutip=6.9083;% for 128x128 mesh

Reference_Sze =[
0.05 0.855 0.955 
0.10 1.499 1.840 
0.15 1.969 2.604 
0.20 2.321 3.261 
0.25 2.596 3.833 
0.30 2.819 4.339 
0.35 3.002 4.790
0.40 3.158 5.196 
0.45 3.291 5.565 
0.50 3.406 5.902 
0.55 3.508 6.212 
0.60 3.598 6.497 
0.65 3.678 6.761 
 0.70 3.750 7.006
0.75 3.816 7.234
0.80 3.875 7.448
0.85 3.929 7.647
0.90 3.979 7.835
0.95 4.025 8.011
1.00 4.067 8.178]';

Reference_Sze_short =[
0.05 0.855 0.955 
0.10 1.499 1.840 
0.15 1.969 2.604 
0.20 2.321 3.261 
0.25 2.596 3.833 
0.30 2.819 4.339 
0.35 3.002 4.790
0.40 3.158 5.196 
0.45 3.291 5.565 
0.50 3.406 5.902 ]';
Reference_Sze_short(1,:) =Reference_Sze_short(1,:)*200;


% 4 x 4 mesh
Mesh_1 =[
    0    10    20    30    40    50    60    70    80    90   100
    0    0.8980    1.5381    1.9762    2.2877    2.5216    2.7054    2.8550    2.9799    3.0866    3.1791
    0   -1.0794   -2.1367   -3.0371   -3.7721   -4.3757   -4.8805   -5.3105   -5.6825   -6.0087   -6.2980];


% 8 x 8 mesh
Mesh_2 =[
    0    10    20    30    40    50    60    70    80    90   100
     0    0.8857    1.5198    1.9847    2.3166    2.5793    2.7798    2.9320    3.0854    3.2002    3.3022
        0   -1.0119   -1.9256   -2.7284   -3.3822   -3.9605   -4.4375   -4.8207   -5.2462   -5.5767   -5.8854
    ];

% 16 x 16 mesh
Mesh_3 =[
    0    10    20    30    40    50    60    70    80    90   100
    0    0.8881    1.5421    2.0132    2.3658    2.6401    2.8572    3.0330    3.1852    3.3112    3.4211
    0   -0.9987   -1.9099   -2.6917   -3.3639   -3.9491   -4.4560   -4.8984   -5.3081   -5.6668   -5.9958
    ];



figure('Position',[100,100,640,600]);

plot(Reference_Sze_short(2,:),Reference_Sze_short(1,:),...
    name_to_style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(-Reference_Sze_short(3,:),Reference_Sze_short(1,:),...
    name_to_style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','Sze','TextBackground','w','color','k','fontname','times','fontsize',36)

plot(Mesh_1(2,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(Mesh_1(3,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,4 x 4','TextBackground','w','color','k','fontname','times','fontsize',36)

plot(Mesh_2(2,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(Mesh_2(3,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,8 x 8','TextBackground','w','color','k','fontname','times','fontsize',36)


plot(Mesh_3(2,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(Mesh_3(3,:),Mesh_2(1,:),...
    name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,16 x 16','TextBackground','w','color','k','fontname','times','fontsize',36)

set(gca,'XLim', [-6.5, 3.5]);
% set(gca,'YLim', [0.001, 1.0]);
xlabel ('Radial deflection')
% ylabel ('Est. True Error of Deflection')
ylabel ('Loading parameter \lambda')
% legend(Description,'Location','Southeast')
set(gca,'linewidth',2);
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);
% set_decades_on_axis (gca)
hold on; grid on; figure (gcf);
% title( ['Fiber reinforced cantilever, iso'])
% saveas(gcf,[mfilename '.png']);
% saveas(gcf,[mfilename '.fig']);
%  saveas(gcf,[mfilename '.eps']);
%  close  all


%
