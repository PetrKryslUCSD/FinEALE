function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.94;% for 128x128 mesh

function s=n2style(d)
switch d
    case 'F-Bar SN Q'
        s='ks-';
    case 'F-Bar SN T'
        s='kv--';
    case 'Q1E4'
        s='kd-';
    case 'T10MS'
        s='rx-';
    otherwise
        s=name_to_style(d);
end
end

Data{end+1}=[
  [  2.0   3.8315
    4.0    6.4359
    8.0    6.7527
   16.0539    6.9073
   32.0    6.9227]'
];'';Description{end+1} ='Q1E4';

Data{end+1}=[
  % H8MSGSO
[2.0  7.5731
4.0   6.0748
8.0 6.9232
16.0   6.8500
32.0   6.9334]'
   ];'';Description{end+1} ='H8MSGSO';

Data{end+1}=[
  % H8MSGSO
[2.0   5.1052
4.0     6.6374
8.0  6.7336
16.0     6.8826
32.0   6.933]'
   ];'';Description{end+1} ='H20R';

% for j=1:(length(Data))
%     data=Data{j};
%     n=44./data(1,:); uz=data(2,:);
%     %     [xestim, beta, c, residual] = richextrapol(se([1,2,3]),h([1,2,3]));
%     %         loglog(n,abs((uz-convutip)/convutip),name_to_style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
%     plot(n,uz/convutip,n2style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
%     annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string',Description{j},'TextBackground','w','color','k','fontname','times','fontsize',36)
% end
%
% % set(gca,'XLim', [100, 1000]);
% % set(gca,'YLim', [0.001, 1.0]);
% xlabel ('Number of elements per side')
% % ylabel ('Est. True Error of Deflection')
% ylabel ('Deflection')
% % legend(Description,'Location','Southeast')
% set(gca,'Position',[0.2 0.17 0.75 0.78]);
% options.FontSize=30;
% set_pub_defaults(gcf,options);
% % set_decades_on_axis (gca)
% hold on; grid on; figure (gcf);
% % title( ['Fiber reinforced cantilever, iso'])
% % saveas(gcf,[mfilename '.png']);
% % saveas(gcf,[mfilename '.fig']);
% %  saveas(gcf,[mfilename '.eps']);
% %  close  all


%
figure('Position',[100,100,640,700]);
for j=1:(length(Data))
    data=Data{j};
    n=data(1,:); uz=data(2,:);
    %     [xestim, beta, c, residual] = richextrapol(se([1,2,3]),h([1,2,3]));
    %             loglog(n,abs((uz-convutip)/convutip),n2style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
                plot(n,uz,name_to_style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string',Description{j},'TextBackground','w','color','k','fontname','times','fontsize',36)
end

set(gca,'XLim', [0, 40]);
% set(gca,'YLim', [0.001, 1.0]);
xlabel ('Number of elements per side')
ylabel ('Est. True Error of Deflection')
% ylabel ('Deflection')
% legend(Description,'Location','Southeast')
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);
% set_decades_on_axis (gca)
hold on; grid on; figure (gcf);
% title( ['Fiber reinforced cantilever, iso'])
% saveas(gcf,[mfilename '.png']);
% saveas(gcf,[mfilename '.fig']);
%  saveas(gcf,[mfilename '.eps']);
%  close  allend
end
