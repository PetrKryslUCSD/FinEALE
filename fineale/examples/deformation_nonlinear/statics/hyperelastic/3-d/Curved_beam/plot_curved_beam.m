Fmag=600;


%%
% Comparison with reference solution.

izzelna_ux =[    0.0317    0.6048
    -0.0628    0.6048
    8.7273   34.7095
    16.6668   85.8664
    23.7557  134.8918
    32.6404  213.7587
    37.8389  264.9157
    41.7141  318.2041
    44.9277  382.1503
    47.9523  446.0965
    50.4098  510.0426
    51.9220  578.2519
    53.9069  657.1188
    53.9069  657.1188];

izzelna_uy =[-0.1574    0.6048
    0.0317    2.7364
    1.1659   73.0772
    2.7727  132.7603
    4.9466  211.6272
    7.1205  284.0995
    8.9164  350.1772
    10.2396  414.1234
    11.3738  480.2011
    13.1697  548.4104
    14.0203  618.7511
    14.6820  657.1188
    14.5875  657.1188];




izzelna_uz =[   -0.5354    2.7364
    -0.2519    2.7364
    -1.1970   51.7618
    -2.2367   90.1295
    -4.3161  124.2341
    -6.8681  183.9172
    -9.1366  235.0741
    -11.8776  288.3626
    -14.5241  354.4403
    -16.9815  411.9919
    -19.3445  475.9380
    -21.7074  550.5419
    -23.9759  618.7511
    -25.3936  661.3819];


figure
plot(izzelna_ux(:,1),izzelna_ux(:,2),'r--','linewidth',2)
hold on
plot(izzelna_uy(:,1),izzelna_uy(:,2),'r--','linewidth',2)
hold on
plot(izzelna_uz(:,1),izzelna_uz(:,2),'r--','linewidth',2)
hold on


lambdas =[   0.1000    0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000    1.0000];

% 8 el lengthwise, 2x2 crossx
u1s =[
    0.5734   -0.8923   11.5785
    2.0308   -3.1887   21.6383
    3.8865   -6.1653   29.6458
    5.7924   -9.2774   35.8078
    7.5820  -12.2489   40.5392
    9.1985  -14.9755   44.2165
   10.6357  -17.4358   47.1241
   11.9071  -19.6432   49.4654
   13.0323  -21.6232   51.3836
   14.0311  -23.4033   52.9805];
tip_displacement=u1s;load_parameters=lambdas;
sty='kx';

plot(tip_displacement(:,1),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)
plot(tip_displacement(:,2),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)
plot(tip_displacement(:,3),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)

% 2 el lengthwise, 1x1 crossx
u1s =[
   2.0914   -0.0900   11.2823
    3.8693   -1.0892   18.3873
    5.5954   -2.4140   24.0027
    7.1956   -3.8596   28.5137
    8.6402   -5.3227   32.1931
    9.9292   -6.7525   35.2418
   11.0758   -8.1254   37.8065
   12.0969   -9.4317   39.9938
   13.0094  -10.6693   41.8823
   13.8284  -11.8395   43.5308];
tip_displacement=u1s;load_parameters=lambdas;
sty='ko';

plot(tip_displacement(:,1),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)
plot(tip_displacement(:,2),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)
plot(tip_displacement(:,3),load_parameters*Fmag,sty,'linewidth',3, 'markersize',8)

% % 4 el lengthwise, 1x1 crossx
% u1s =[
%    3.4876   -1.1304   17.3536
%     5.7996   -3.5363   25.7272
%     7.7522   -6.1938   31.8493
%     9.4320   -8.7919   36.4978
%    10.8850  -11.2151   40.1161
%    12.1488  -13.4310   42.9940
%    13.2549  -15.4419   45.3275
%    14.2295  -17.2632   47.2523
%    15.0937  -18.9144   48.8643
%    15.8649  -20.4153   50.2327];
% tip_displacement=u1s;load_parameters=lambdas;
% sty='kd';
%
% plot(tip_displacement(:,1),load_parameters*Fmag,sty,'linewidth',3)
% plot(tip_displacement(:,2),load_parameters*Fmag,sty,'linewidth',3)
% plot(tip_displacement(:,3),load_parameters*Fmag,sty,'linewidth',3)


% title('Curved cantilever, L-to-R uy,ux,uz','fontname','Times','fontsize',16)

set(gca,'XLim', [-30, 60]);
set(gca,'YLim', [0,600]);
labels(' Tip displacement [in]','Load [lb]')
grid on
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);

