font='Times';
Legends = {};
Data = {};
Style = {};
set_graphics_defaults

% loadv=[0;0;p];dir=3;uzex=5.426e-3;
% loadv=[0;p;0];dir=2;
%uzex=15.41;
ref = [22, 42.5, 59, 70];
    
% abaqus
Data{end+1}=[
    [2, 8, 16];
  ref',ref',ref'
]; Style{end+1}='k-'; Legends{end+1} ='Ref';

Data{end+1}=[
    [4, 8, 16];
  [66.5,66.5,69];;[66.5,66.5,69];[66.5,66.5,69];[66.5,66.5,69]
]; Style{end+1}='md-'; Legends{end+1} ='CayMah';

Data{end+1}=[
    [4, 8, 16];
  [66.5,72.3,71];;[66.5,72.3,71];[66.5,72.3,71];[66.5,72.3,71]
]; Style{end+1}='bv-'; Legends{end+1} ='Masud';

% Unstructured,, tetrahedral  to hexahedral,, grid stabilization fraction 0.75
Data{end+1}=[
    [2,4,8,16];
    [0.25, 0.45, 0.66, 0.78
    0.225, 0.41, 0.576,  0.705
    00.22, 0.43, 0.58, 0.69
    0.22, 0.43, 0.59, 0.695
   ]'*100;
]; Legends{end+1} ='Q1SP'; Style{end+1}='r-x'; 


Data{end+1}=[
    [4 8];
  [0.2040  0.2150]*100;;[  0.3875   0.4127 ]*100;[ 0.5281  0.5663  ]*100;[  0.6247 0.6706]*100
]; Style{end+1}='r^-'; Legends{end+1} ='T10MS';

% % Regular grid, stabilization fraction 0.75
% Data{end+1}=[
%     [2,4,6,8,10];
%     [0.248946188027685   0.497699002324237   0.726791486473527   0.915696437451734
%     0.218230580047313   0.422360867266289   0.604924858125240   0.772583918821509
%     0.222263737332980   0.431419027723265   0.603834120417362   0.746580090778514
%     0.219616809652744   0.425970693078486   0.592463823354558   0.720334557027934
%     0.220402769607209   0.426661435811957   0.590155104215669   0.709871451321210
%    ]'*100;
% ]; Legends{end+1} ='H8MSGS'; Style{end+1}=name_to_style(Legends{end}); 
% 
% % Unstructured,, tetrahedral  to hexahedral,, grid stabilization fraction 0.75
% Data{end+1}=[
%     [4,6,8,12];
%     [0.213017898365441   0.409365294966901   0.556585311572643   0.658049884776807
%     0.218895961938956   0.425128216989716   0.587486541943935   0.700017392690646
%     0.221716723715336   0.433195194840638   0.601661607915919   0.717842724192794
%     2.203998840491740e-01     4.271931761111313e-01     5.900305960529763e-01     7.004029071585533e-01
%    ]'*100;
% ]; Legends{end+1} ='H8MSGS(U)'; Style{end+1}='k-v'; 

% Unstructured,, tetrahedral  to hexahedral,, grid stabilization fraction
% optimal
% Data{end+1}=[
%     [4,8,12];
%     [0.2128    0.4087    0.5586    0.6652
%     0.2217    0.4300    0.5950    0.7094
%     0.2200    0.4261    0.5878    0.6997
%    ]'*100;
% ]; Legends{end+1} ='H8MSGSO(U)'; Style{end+1}='k-v'; 

% Structured grid stabilization fraction
% optimal
Data{end+1}=[
    [4,8,12];
    [ 0.2059    0.3951    0.5444    0.6517
     0.2213    0.4277    0.5903    0.7042
     0.2180    0.4211    0.5809    0.6914
    ]'*100;
]; Legends{end+1} ='H8MSGSO'; Style{end+1}='ko-'; 



% Data{end+1}=[
%     [ 2, 4, 6];
%    [0.1895, 0.211, 0.21572]*100;
%   [0.3365, 0.41908,0.42935]*100;
%   [0.4528,0.57095, 0.5729]*100;
%   [0.5272, 0.67165,0.697]*100;
% ]; Style{end+1}='mh-'; Legends{end+1} ='NICE-H27';


% Data{end+1}=[
%     [2,4,6];
%    [0.203, 0.21348, 0.21679]*100;
%   [0.354, 0.39687,0.4101]*100;
%   [0.442442,0.5217, 0.55423]*100;
%   [0.49656, 0.6165,0.645]*100;
% ]; Style{end+1}='ks-'; Legends{end+1} ='H27';

for k=2:5
    for j=1:length(Data)
        h=Data{j}(1,:);
        data=Data{j}(k,:);
        plot(h,data,Style{j},'lineWIDTH', 3,'markersize',8); hold on;
    end
    %     loglog(h,abs(utip-uzex)/uzex,Style{j},'lineWIDTH', 2); hold on;
end

set(gca,'XLim', [2,16]);
% set(gca,'YLim', [0.4, 1.3]);
% set(gca,'YLim', [0.9, 1.05]);
xlabel ('Number of elements per side')
ylabel ('Compression [%]')
legend(Legends,'Location','Southeast')
set(gca,'Position',[0.15 0.15 0.5 0.5]);
% set(gca,'XTickLabels', {'u_x','','u_y','','u_z'});
% set_presentation_defaults
options.FontSize=30;
set_pub_defaults(gcf,options);
set(gca,'lineWIDTH', 2)
hold on; grid on; figure (gcf);
% saveas(gcf,[mfilename '.png']);
 printeps(gcf,[mfilename '.eps']);
%  close  all