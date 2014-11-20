function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.932;
convutip=6.9083;% for 128x128 mesh

function s=n2style(d)
switch d
    case 'SHELL'
        s='k^-.';
    case 'Q1SP'
        s='kd-';
    case 'Q1M/E12'
        s='k+--';
    case '1 el'
        s='ko--';
    case '2 el'
        s='ko:';
    case '4 el'
        s='ko-';
    otherwise
        s=name_to_style(d);
end
end

Data{end+1}=[
  8  16  32
  -12.7 -10.7 -14.3

];'';Description{end+1} ='Q1SP';
Data{end+1}=[
  8  16  32
  -2.7 -8.1 -14.

];'';Description{end+1} ='Q1M/E12';

Data{end+1}=[
  8  16  32
  -10 -14 -15.5

];'';Description{end+1} ='SHELL';

%    0    2.1250    4.2500    6.3750    8.5000
%   0   -2.9497   -7.2015  -11.3577  -15.2371
%  0   -2.6041   -6.0674   -9.5259  -12.2170
Data{end+1}=[
8 12 16 
  -12.2497 -13.9851 -14.78
];'';Description{end+1} ='4 el';
Data{end+1}=[
 8 12 16 
 -12.25 -14.4375 -15.3068
];'';Description{end+1} ='2 el';

% 1 el / thickness
Data{end+1}=[
  8 12 16.0000 
  -13.27 -16.486 -17.587
];'';Description{end+1} ='1 el';

figure('Position',[100,100,640,700]);
for j=1:(length(Data))
    data=Data{j};
    n=data(1,:); uz=-data(2,:);
    %     [xestim, beta, c, residual] = richextrapol(se([1,2,3]),h([1,2,3]));
    %         loglog(n,abs((uz-convutip)/convutip),name_to_style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
    plot(n,uz,n2style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string',Description{j},'TextBackground','w','color','k','fontname','times','fontsize',36)
end
%
set(gca,'XLim', [8, 32]);
set(gca,'YLim', [0.0, 18]);
xlabel ('Number of elements circumferentially')
% ylabel ('Est. True Error of Deflection')
ylabel ('Deflection at A')
% legend(Description,'Location','Southeast')
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);
set(gca,'lineWIDTH', 2)
hold on; grid on; figure (gcf);
% % title( ['Fiber reinforced cantilever, iso'])
% % saveas(gcf,[mfilename '.png']);
% % saveas(gcf,[mfilename '.fig']);
% %  saveas(gcf,[mfilename '.eps']);
% %  close  all


%
% figure('Position',[100,100,640,700]);
% for j=1:(length(Data))
%     data=Data{j};
%     n=44./data(1,:); uz=data(2,:);
%     %     [xestim, beta, c, residual] = richextrapol(se([1,2,3]),h([1,2,3]));
%             loglog(n,abs((uz-convutip)/convutip),n2style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
%             %     plot(n,uz,name_to_style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
%     annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string',Description{j},'TextBackground','w','color','k','fontname','times','fontsize',36)
% end
%
% % set(gca,'XLim', [100, 1000]);
% % set(gca,'YLim', [0.001, 1.0]);
% xlabel ('Number of elements per side')
% ylabel ('Est. True Error of Deflection')
% % ylabel ('Deflection')
% % legend(Description,'Location','Southeast')
% set(gca,'Position',[0.2 0.17 0.75 0.78]);
% options.FontSize=30;
% set_pub_defaults(gcf,options);
% set_decades_on_axis (gca)
% hold on; grid on; figure (gcf);
% % title( ['Fiber reinforced cantilever, iso'])
% % saveas(gcf,[mfilename '.png']);
% % saveas(gcf,[mfilename '.fig']);
% %  saveas(gcf,[mfilename '.eps']);
% %  close  allend
% end

end