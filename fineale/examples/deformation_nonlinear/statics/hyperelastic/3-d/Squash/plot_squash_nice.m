font='Times';
Legends = {};
Data = {};
Style = {};
set_graphics_defaults

% loadv=[0;0;p];dir=3;uzex=5.426e-3;
% loadv=[0;p;0];dir=2;
%uzex=15.41;
ref = [22, 42.5, 59, 70];
h=[4, 6, 8];
    
% abaqus
Data{end+1}=[
    [2, 6, 8];
  ref',ref',ref'
]; Style{end+1}='k--'; Legends{end+1} ='Ref';

Data{end+1}=[
    [4, 6, 8];
   [0.2003, 0.2035, 0.2076]*100;
  [0.4021, 0.4176,0.4308]*100;
  [0.5641,0.5943, 0.6120]*100;
  [0.6683, 0.7177,0.7337]*100;
]; Style{end+1}='m<-.'; Legends{end+1} ='NICE-T4';

Data{end+1}=[
    [ 2, 4, 6];
   [0.1895, 0.211, 0.21572]*100;
  [0.3365, 0.41908,0.42935]*100;
  [0.4528,0.57095, 0.5729]*100;
  [0.5272, 0.67165,0.697]*100;
]; Style{end+1}='mh-'; Legends{end+1} ='NICE-H27';


Data{end+1}=[
    [2,4,6];
   [0.203, 0.21348, 0.21679]*100;
  [0.354, 0.39687,0.4101]*100;
  [0.442442,0.5217, 0.55423]*100;
  [0.49656, 0.6165,0.645]*100;
]; Style{end+1}='ks-'; Legends{end+1} ='H27';

for k=2:5
    for j=1:length(Data)
        h=Data{j}(1,:);
        data=Data{j}(k,:);
        plot(h,data,Style{j},'lineWIDTH', 2); hold on;
    end
    %     loglog(h,abs(utip-uzex)/uzex,Style{j},'lineWIDTH', 2); hold on;
end

% set(gca,'XLim', [100, 1000]);
% set(gca,'YLim', [0.4, 1.3]);
% set(gca,'YLim', [0.9, 1.05]);
xlabel ('Number of elements per side')
ylabel ('Compression [%]')
legend(Legends,'Location','Southeast')
set(gca,'Position',[0.15 0.15 0.5 0.5]);
% set(gca,'XTickLabels', {'u_x','','u_y','','u_z'});
set_presentation_defaults
hold on; grid on; figure (gcf);
% saveas(gcf,[mfilename '.png']);
 printeps(gcf,[mfilename '.eps']);
%  close  all