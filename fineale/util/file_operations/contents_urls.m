
function contents_urls()
    sep=filesep;
    dl=dir(pwd);
    fid=fopen('url.asv','w');
    for i=1:length(dl)
        name = strrep([strrep(pwd,[fineale_path filesep],'') filesep dl(i).name],'\','/')
        
        if (dl(i).isdir)
            if      (~strcmp(dl(i).name,'.')) && ...
                    (~strcmp(dl(i).name,'..'))
                % <a href="matlab: doc 'fahy/Contents.m'">fahy</a> - one-dimensional examples
                [pathstr, name1, ext] = fileparts(name);
                fprintf(fid,['%% <a href="matlab: doc ''' name '/Contents.m''">' name1 '</a> - \n']);
            end
        else
            if      (~strcmp(dl(i).name,'url.asv'))&& ...
                    (~strcmp(dl(i).name,'Contents.m'))
                %   <a href="matlab: edit 'underintegr2'">underintegr2</a>     - Example with underintegration of Q4
                [pathstr, name1, ext] = fileparts(name);
                fprintf(fid,['%% <a href="matlab: edit ''' name '''">' name1 '</a> - \n']);
            end
        end
    end
    fclose(fid);
end
% function contents_urls()
%     sep=filesep;
%     dl=dir(pwd);
%     for i=1:length(dl)
% %         dl(i)
% disp(dl(i).name)
%         if (dl(i).isdir)
%             if      (~strcmp(dl(i).name,'.')) && ...
%                     (~strcmp(dl(i).name,'..'))
%                 run_in_subdirs([d sep dl(i).name], Suffixes, apply);
%             end
%         else
%             [pathstr, name, ext] = fileparts(dl(i).name) ;
%             for j=1: length(Suffixes)
%                 if ((strcmp(ext,Suffixes{j})))
%                     apply([d sep dl(i).name]);
%                 end
%             end
%         end
%     end
% end
