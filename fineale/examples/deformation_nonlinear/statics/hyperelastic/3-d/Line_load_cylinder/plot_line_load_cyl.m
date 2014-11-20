function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.932;
convutip=6.9083;% for 128x128 mesh

function s=n2style(d)
switch d
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
  -13.75 -15.6 -16.4

];'';Description{end+1} ='Q1SP';
Data{end+1}=[
  8  16  32
  -12.3 -16.05 -16.4

];'';Description{end+1} ='Q1M/E12';


%  eight elements circumferentially
% Data{end+1}=[
%    0        1875        3750        5625        7500
% 0   -2.6113   -6.0320   -9.7988  -13.1284
% ];'';Description{end+1} ='8';
% %  12 elements circumferentially
% Data{end+1}=[
%    0        1875        3750        5625        7500
% 0   -2.7860   -6.6117  -11.0482  -14.9695
% ];'';Description{end+1} ='12';
% 16 elements
% 0   -2.8309   -6.7804  -11.4381  -15.5720
% 20 elements
% 0   -2.8442   -6.8381  -11.5858  -15.8162
% 32 elements
% 0   -2.9282   -7.0939  -12.0450  -16.4135
Data{end+1}=[
  8 12 16 
   -13.1845 -14.9926  -15.7022
];'';Description{end+1} ='4 el';
Data{end+1}=[
 8  12 16 
 -12.8708 -14.7057 -15.3035
];'';Description{end+1} ='2 el';
Data{end+1}=[
  8 12 16 
  -12.4609 -13.7576 -14.9255
];'';Description{end+1}  ='1 el';

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
set(gca,'YLim', [10, 17]);
xlabel ('Number of elements circumferentially')
% ylabel ('Est. True Error of Deflection')
ylabel ('Deflection at A')
% legend(Description,'Location','Southeast')
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);
set(gca,'lineWIDTH', 2)
hold on; grid on; figure (gcf);


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