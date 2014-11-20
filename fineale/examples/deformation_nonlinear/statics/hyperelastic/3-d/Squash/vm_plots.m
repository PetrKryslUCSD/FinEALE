set_graphics_defaults
set_label_defaults

for i=1:10:length(us) 
    u1=us{i};
    vm_plot
      print(gcf, '-r50', [mfilename n2s0p(i) '.png'], '-dpng');
end
    