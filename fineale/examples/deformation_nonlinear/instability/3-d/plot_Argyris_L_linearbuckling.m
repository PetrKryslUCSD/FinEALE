font='Times';
Legends = {};
Data = {};
Style = {};
set_graphics_defaults

%%
% References
% J.H. Argyris, P.C. Dunne and D.W. Scharpf, On large
% displacement-small strain analysis of structures with
% rotation degree of freedom. Comput. Methods Appl. Mech. Engrg. 14 (1978) 401-451; 15 (1978) 99-135.
% J.H. Argyris, 0. Hilpert, GA. Malejannakis and D.W. Scharpf, On the
% geometrical stiffness of a beam in space-A consistent V.W.
% Approach, Comput. Methods Appl. Mech. Engrg. 20 (1979) 105-131.
%
% The reference values of the buckling parameters are 68040 (for the load
% putting the clamped leg into compression) and -108800 (for the load
% putting the clamped leg in tension).
ref = [-68040,108800];
    
% abaqus
Data{end+1}=[
    [1, 9];
  ref',ref'
]; Style{end+1}='k-'; Legends{end+1} ='Ref';

Data{end+1}=[
1                     2                     3                     4                     5                     6                     7                     8                     9
85592.6487322732      93467.1970399187       112449.99962936      113161.797625028       116471.65818484      116688.319467014      118561.846864718      119149.566793323      120766.716399375
-126531.600898622     -77101.4867171661     -74108.1071911146     -72136.5140903131     -71390.4539144491      -71467.639391695      -71196.738891624     -71525.7017658502     -71562.0771230718
];Style{end+1}='ko-';Legends{end+1} ='H8MSGSO';
Data{end+1}=[
1                     2                     3                     4                     5                     6                     7                     8                     9
132180.592735386      125178.632656284        122392.3220599       121080.06583636       120393.75719104      120016.522023864      119742.906689001      119603.359320247      119460.314328464
-76707.5021648862      -72738.703086779     -71781.7821786533     -71287.5739389513     -71035.2393315045     -70873.8587820839     -70745.4865180825     -70693.8955565148     -70607.9683870395
];Style{end+1}='kp-';Legends{end+1} ='H64';
Data{end+1}=[
1                     2                     3                     4                     5                     6                     7                     8                     9
158733.875700293      137225.951114168      132016.476145931       129048.69015239      127350.185598022      126151.256224438      125280.788304616      124638.413613235      124136.045442173
-108547.001358987      -84482.103639827      -78466.196414641     -75891.9824777608     -74436.5204446176     -73555.9706831883     -72950.5006361363     -72513.0112280041     -72195.1519533478
];Style{end+1}='cd-';Legends{end+1} ='H27';

 for j=1:length(Data)
        h=Data{j}(1,:);
        data=Data{j}(2,:);
        plot(h,data,Style{j},'lineWIDTH', 2,'markersize',8); hold on;
    end

 for j=1:length(Data)
        h=Data{j}(1,:);
        data=Data{j}(3,:);
        plot(h,data,Style{j},'lineWIDTH', 2,'markersize',8); hold on;
    end

set(gca,'XLim', [1,9]);
% set(gca,'YLim', [0.4, 1.3]);
% set(gca,'YLim', [0.9, 1.05]);
xlabel ('Number of elements per side [ND]')
ylabel ('Buckling load [N]')
legend(Legends,'Location','Southeast')
set(gca,'Position',[0.15 0.15 0.5 0.5]);
% set(gca,'XTickLabels', {'u_x','','u_y','','u_z'});
options.FontSize=30;
set_pub_defaults(gcf,options);
hold on; grid on; figure (gcf);
% saveas(gcf,[mfilename '.png']);
%  printeps(gcf,[mfilename '.eps']);
%  close  all