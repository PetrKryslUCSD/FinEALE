function plot_nonlinear_twisted_beam
font='Times';
Legends = {};
Data = {};
Description = {};

refx=[0,0
    0.1078    1.8251
    0.2994    3.5361
    0.6707    5.4753
    1.1018    7.3004
    1.6766    9.4677
    2.2275   11.2928
    2.6587   13.3460
    3.1138   14.9430
    3.4850   16.6540
    4.1198   19.7338
    4.9820   25.2091
    5.7246   31.1407
    6.4072   37.7567
    6.7904   43.0038
    6.9341   45.3992
    7.3511   54.2966
    7.6048   60.1141];
refy=[0 0
    0.4198    1.5849
    0.9355    3.1698
    1.8591    5.5472
    2.8186    8.3774
    3.6582   10.8679
    4.2819   13.4717
    4.5817   14.6038
    4.8096   15.9623
    5.0375   16.9811
    5.2894   18.5660
    5.6972   21.7358
    6.0810   25.1321
    6.4048   28.9811
    6.6927   33.0566
    6.9685   37.8113
    7.2204   43.0189
    7.4123   48.3396
    7.6042   54.2264
    7.7961   60.0000];
refz=[0 0
    1.3154    3.1638
    1.9851    4.5198
    2.6069    6.3277
    3.1689    7.6836
    3.7668    9.9435
    4.0897   11.7514
    4.3408   13.5593
    4.5082   15.1412
    4.6398   16.9492
    4.7115   20.0000
    4.7593   25.3107
    4.6637   31.0734
    4.5441   37.8531
    4.3886   45.5367
    4.2332   54.4633
    4.1256   59.8870];

function s=n2style(d)
switch d
    case '4x2x2'
        s='kd';
    case '8x4x4'
        s='kv';
    otherwise
        s=name_to_style(d);
end
end

lam=[     6    12    18    24    30    36    42    48    54    60];
Mesh_4_2_2=[
   -0.7646    2.1394   -2.7885
   -2.2681    4.0478   -4.0789
   -3.4292    5.1616   -4.3659
   -4.2318    5.8152   -4.3635
   -4.8085    6.2355   -4.2848
   -5.2426    6.5275   -4.1902
   -5.5823    6.7432   -4.0982
   -5.8565    6.9101   -4.0141
   -6.0835    7.0442   -3.9388
   -6.2753    7.1552   -3.8720
];

Mesh_8_4_4=[
    -0.7438    2.0846   -2.7152
   -2.2939    3.9104   -4.1792
   -3.6637    5.1009   -4.6642
   -4.6920    5.8587   -4.7694
   -5.4640    6.3744   -4.7386
   -6.0585    6.7491   -4.6574
   -6.5280    7.0355   -4.5584
   -6.9072    7.2633   -4.4554
   -7.2195    7.4499   -4.3544
   -7.4809    7.6066   -4.2581   ];

figure('Position',[100,100,640,600]);

% plot(Reference_Sze_short(2,:),Reference_Sze_short(1,:),...
%     name_to_style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
% plot(-Reference_Sze_short(3,:),Reference_Sze_short(1,:),...
%     name_to_style('reference'),'lineWIDTH', 3, 'markersize',10); hold on;
%     annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','Sze','TextBackground','w','color','k','fontname','times','fontsize',36)

plot(refx(:,1),refx(:,2),...
    n2style('reference'),'lineWIDTH', 3, 'markersize',5); hold on;
plot(refy(:,1),refy(:,2),...
    n2style('reference'),'lineWIDTH', 3, 'markersize',5); hold on;
plot(refz(:,1),refz(:,2),...
    n2style('reference'),'lineWIDTH', 3, 'markersize',5); hold on;

plot(-Mesh_4_2_2(:,1),lam,...
    n2style('4x2x2'),'lineWIDTH', 3, 'markersize',8); hold on;
plot(Mesh_4_2_2(:,2),lam,...
    n2style('4x2x2'),'lineWIDTH', 3, 'markersize',8); hold on;
plot(-Mesh_4_2_2(:,3),lam,...
    n2style('4x2x2'),'lineWIDTH', 3, 'markersize',8); hold on;
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,4 x 4','TextBackground','w','color','k','fontname','times','fontsize',36)

plot(-Mesh_8_4_4(:,1),lam,...
    n2style('8x4x4'),'lineWIDTH', 3, 'markersize',8); hold on;
plot(Mesh_8_4_4(:,2),lam,...
    n2style('8x4x4'),'lineWIDTH', 3, 'markersize',8); hold on;
plot(-Mesh_8_4_4(:,3),lam,...
    n2style('8x4x4'),'lineWIDTH', 3, 'markersize',8); hold on;
annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string','H8MSGSO,4 x 4','TextBackground','w','color','k','fontname','times','fontsize',36)

set(gca,'XLim', [0, 8]);
set(gca,'YLim', [0.00, 60]);
xlabel ('Tip deflection')
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
