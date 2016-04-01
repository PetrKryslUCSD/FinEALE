function plot_two_layer
font='Times';
Legends = {};
Data = {};
Description = {};


function s=n2style(d)
switch d
    case '32x12x4'
        s='ks-';
    case '16x6x2'
        s='kd--';
    case '8x3x1'
        s='ko-.';
    case '8x3x1_H64'
        s='rh-.';
    otherwise
        s=name_to_style(d);
end
end

% 4 x 4 x 2 mesh
Mesh_32x12x4 =[
0.0029   -0.0001    0.0077
    0.0057    0.0006    0.0127
    0.0084    0.0014    0.0162
    0.0110    0.0022    0.0189
    0.0136    0.0030    0.0210
    0.0163    0.0038    0.0227
    0.0188    0.0045    0.0242
    0.0214    0.0052    0.0255
    0.0240    0.0059    0.0266
    0.0265    0.0065    0.0276
    ];
% % 8 x 8 x 3 mesh
Mesh_16x6x2 =[
  0.0029   -0.0000    0.0082
    0.0057    0.0007    0.0134
    0.0084    0.0016    0.0170
    0.0110    0.0025    0.0197
    0.0137    0.0033    0.0219
    0.0163    0.0041    0.0237
    0.0189    0.0049    0.0251
    0.0214    0.0056    0.0264
    0.0240    0.0063    0.0276
    0.0265    0.0070    0.0286  ];
%
% % 8 x 8 x 4 mesh
Mesh_8x3x1 =[
0.0030    0.0002    0.0109
    0.0058    0.0014    0.0170
    0.0086    0.0027    0.0209
    0.0113    0.0039    0.0238
    0.0139    0.0050    0.0260
    0.0166    0.0059    0.0277
    0.0192    0.0069    0.0292
    0.0218    0.0077    0.0304
    0.0244    0.0086    0.0315
    0.0270    0.0093    0.0325];

Mesh_8x3x1_H64=[    0.0029   -0.0001    0.0075
    0.0056    0.0006    0.0124
    0.0083    0.0013    0.0159
    0.0110    0.0021    0.0185
    0.0136    0.0029    0.0207
    0.0162    0.0036    0.0224
    0.0188    0.0044    0.0239
    0.0214    0.0050    0.0251
    0.0239    0.0057    0.0262
    0.0265    0.0063    0.0272];

figure('Position',[100,100,640,600]);

% plot(Reference_Sze_short(2,:),Reference_Sze_short(1,:),...
%     n2style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
% plot(-Reference_Sze_short(3,:),Reference_Sze_short(1,:),...
%     n2style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
%     annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','Sze','TextBackground','w','color','k','fontname','times','fontsize',36)

ddata =Mesh_32x12x4;
for k=1:3
    plot([0;ddata(:,k)],0:1:10',...
        n2style('32x12x4'),'lineWIDTH', 3, 'markersize',10); hold on;
end
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','32x12x4','TextBackground','w','color','k','fontname','times','fontsize',36)

ddata =Mesh_16x6x2;
for k=1:3
    plot([0;ddata(:,k)],0:1:10',...
        n2style('16x6x2'),'lineWIDTH', 3, 'markersize',10); hold on;
end
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','16x6x2','TextBackground','w','color','k','fontname','times','fontsize',36)

ddata =Mesh_8x3x1;
for k=1:3
    plot([0;ddata(:,k)],0:1:10',...
        n2style('8x3x1'),'lineWIDTH', 3, 'markersize',10); hold on;
end
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','8x3x1','TextBackground','w','color','k','fontname','times','fontsize',36)

ddata =Mesh_8x3x1_H64;
for k=1:3
    plot([0;ddata(:,k)],0:1:10',...
        n2style('8x3x1_H64'),'lineWIDTH', 3, 'markersize',10); hold on;
end
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','8x3x1_H64','TextBackground','w','color','k','fontname','times','fontsize',36)

% set(gca,'XLim', [-6.5, 3.5]);
% set(gca,'YLim', [0.001, 1.0]);
xlabel ('Deflection of point A')
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

