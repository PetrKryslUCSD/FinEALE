function plot_Cook_3D_strain_abaqus
font='Times';
Legends = {};
Data = {};
Description = {};
convutip=6.932;
convutip=6.9083;% for 128x128 mesh

function s=n2style(d)
switch d
    case 'F-Bar SN Q'
        s='ks-';
    case 'F-Bar SN T'
        s='kv--';
    case 'Q1E4'
        s='kd-';
    otherwise
        s=name_to_style(d);
end
end

Data{end+1}=[
   22.0000   11.0000    2.7500    1.3750 
   [5.332, 6.096, 6.6, 6.661]*1.0392
];'';Description{end+1} ='Q1E4';
Data{end+1}=[
   22.0000   11.0000    5.5000    2.7500    1.3750
   4.847, 6.173, 6.673, 6.833, 6.891
];Description{end+1} ='F-Bar SN Q';
Data{end+1}=[
   22.0000   11.0000    5.5000    2.7500    1.3750
   4.0605, 5.7504, 6.5034, 6.7704, 6.867
];Description{end+1} ='F-Bar SN T';

% Data{end+1}=[
%    5.500000000000000   7.333333333333333  11.000000000000000  22.000000000000000
%    6.895646945396894   6.852368871271234   6.642006177730575   4.883541280078670
%    ];'';Description{end+1} ='H8MSGSO';% for thickness 40.0
% Data{end+1}=[
%     3.1429    3.6667    4.4000    5.5000    7.3333   11.0000   22.0000
%       6.9121    6.9006    6.8784    6.8292    6.6944    6.1879    4.3408
%    ];'';Description{end+1} ='H8MSGSO';% for thickness 20.0
% Data{end+1}=[
%   2.7500    3.1429    3.6667    4.4000    5.5000    7.3333   11.0000   22.0000
%     6.8906    6.8728    6.8429    6.7885    6.6780    6.4214    5.9233    5.0134
%    ];'';Description{end+1} ='H8MSGSO';% for thickness 10.0
Data{end+1}=[  2.7500    3.1429    3.6667    4.4000    5.5000    7.3333   11.0000   22.0000
    6.8627    6.8502    6.8325    6.8056    6.7605    6.6711    6.4273    5.4664
   
   ];'';Description{end+1} ='H8MSGSO';% for thickness 16/n




% Data{end+1}=[
%    22.000000000000000  11.000000000000000   7.333333333333333   5.500000000000000
%    5.752894698151852   6.394590517131232   6.608753451671958   6.705988200436940
% ];'';Description{end+1} ='H20R';

figure('Position',[100,100,640,700]);
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
    n=44./data(1,:); uz=data(2,:);
    %     [xestim, beta, c, residual] = richextrapol(se([1,2,3]),h([1,2,3]));
            loglog(n,abs((uz-convutip)/convutip),n2style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
            %     plot(n,uz,name_to_style(Description{j}),'lineWIDTH', 3, 'markersize',10); hold on;
    annotation('textarrow', [0.5,  0.75], [0.5,  0.75], 'string',Description{j},'TextBackground','w','color','k','fontname','times','fontsize',36)
end

% set(gca,'XLim', [100, 1000]);
% set(gca,'YLim', [0.001, 1.0]);
xlabel ('Number of elements per side')
ylabel ('Est. True Error of Deflection')
% ylabel ('Deflection')
% legend(Description,'Location','Southeast')
set(gca,'Position',[0.2 0.17 0.75 0.78]);
options.FontSize=30;
set_pub_defaults(gcf,options);
set_decades_on_axis (gca)
hold on; grid on; figure (gcf);
% title( ['Fiber reinforced cantilever, iso'])
% saveas(gcf,[mfilename '.png']);
% saveas(gcf,[mfilename '.fig']);
%  saveas(gcf,[mfilename '.eps']);
%  close  allend
end
