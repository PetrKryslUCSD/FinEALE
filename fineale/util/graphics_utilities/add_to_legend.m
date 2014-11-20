% Add a string to the legend
% 
% function handle= add_to_legend(Label)
%     
% If a line with the same label already exists nothing happens; otherwise
% the new label is added to the legend.
% 
function handle= add_to_legend(Label)
    handle=legend('hide');
    s=get(handle,'String');
    for j=1:length(s)
        if (strcmp(s{j},Label))
            handle=legend('show');
            return;;
        end
    end
    s{end+1}=  Label;
    handle =legend(s,'Location','Best');;
    position =get( handle, 'position');
    position (3)= 1.2*position(3);
    set (handle, 'position', position);;
    %     set (handle, 'location', 'best');;
end
