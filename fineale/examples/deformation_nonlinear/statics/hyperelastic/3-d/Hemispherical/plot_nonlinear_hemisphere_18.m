function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.932;
convutip=6.9083;% for 128x128 mesh


function s=n2style(d)
switch d
    case 'Mesh_1'
        s='ks';
    case 'Mesh_2'
        s='kv';
    case 'Mesh_3'
        s='kd';
    otherwise
        s=name_to_style(d);
end
end

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
    0 0 0 
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


% 4 x 4 x 2 mesh
Mesh_1 =[
    0    10    20    30    40    50    60    70    80    90   100
    0    0.9716    1.6326    2.0693    2.3736    2.5988    2.7735    2.9140    3.0302    3.1284    3.2129
0   -1.1779   -2.2898   -3.2058   -3.9411   -4.5394   -5.0362   -5.4563   -5.8174   -6.1318   -6.4088
];
% % 8 x 8 x 3 mesh
Mesh_2 =[
    0    10    20    30    40    50    60    70    80    90   100
    0    0.9064    1.5629    2.0262    2.3668    2.6279    2.8350    3.0037    3.1438    3.2622    3.3637
0   -1.0352   -1.9818   -2.7800   -3.4521   -4.0265   -4.5245   -4.9612   -5.3480   -5.6932   -6.0034
    ];
%
% % 8 x 8 x 4 mesh
Mesh_3 =[
    0    10    20    30    40    50    60    70    80    90   100
 0    0.8914    1.5427    2.0057    2.3473    2.6100    2.8188    2.9891    3.1309    3.2509    3.3538
    0   -1.0151   -1.9481   -2.7389   -3.4073   -3.9798   -4.4771   -4.9141   -5.3017   -5.6482   -5.9599
];



figure('Position',[100,100,640,600]);

plot(Reference_Sze_short(2,:),Reference_Sze_short(1,:),...
    n2style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(-Reference_Sze_short(3,:),Reference_Sze_short(1,:),...
    n2style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','Sze','TextBackground','w','color','k','fontname','times','fontsize',36)

plot(Mesh_1(2,:),Mesh_1(1,:),...
    n2style('Mesh_1'),'lineWIDTH', 3, 'markersize',10); hold on;
plot(Mesh_1(3,:),Mesh_1(1,:),...
    n2style('Mesh_1'),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,4 x 4','TextBackground','w','color','k','fontname','times','fontsize',36)

    plot(Mesh_2(2,:),Mesh_2(1,:),...
        n2style('Mesh_2'),'lineWIDTH', 3, 'markersize',10); hold on;
    plot(Mesh_2(3,:),Mesh_2(1,:),...
        n2style('Mesh_2'),'lineWIDTH', 3, 'markersize',10); hold on;
        annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,8 x 8','TextBackground','w','color','k','fontname','times','fontsize',36)
    %
    %         plot(Mesh_3(2,:),Mesh_3(1,:),...
    %         n2style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    %     plot(Mesh_3(3,:),Mesh_3(1,:),...
    %         n2style('Mesh_3'),'lineWIDTH', 3, 'markersize',10); hold on;
    %         annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,8 x 8','TextBackground','w','color','k','fontname','times','fontsize',36)
    %
    %
    % plot(Mesh_3(2,:),Mesh_2(1,:),...
    %     name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    % plot(Mesh_3(3,:),Mesh_2(1,:),...
    %     name_to_style('H8MSGSO'),'lineWIDTH', 3, 'markersize',10); hold on;
    %     annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,16 x 16','TextBackground','w','color','k','fontname','times','fontsize',36)

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

end
%
