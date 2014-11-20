function sniff_out(Suffixes, apply)
% Search folders and apply some function to matching files.
%
% function sniff_out(Suffixes, apply)
%
% Search the current directory and its subdirectories and execute the apply
% function on the files. 
%
% Suffixes = cell array of suffixes of files to remove (default is 
%           {'.asv','.log','.bak','.m~1~','.m~2~','.m~3~'})
% apply  = handle to a function that takes the name of a file
%
%
% See also: clean_sweep
%
    if  ~exist('Suffixes')
        Suffixes ={'.asv','.log','.bak','.m~1~','.m~2~','.m~3~'};
    end
    function  echo_found(x)
        disp( ['Found ' x]);
    end
    if  ~exist('apply')
        apply =@echo_found;
    end
    run_in_subdirs(pwd, Suffixes, apply)
    return;
end

function run_in_subdirs(d, Suffixes, apply)
    sep=filesep;
    dl=dir(d);
    for i=1:length(dl)
%         dl(i)
        if (dl(i).isdir)
            if      (~strcmp(dl(i).name,'.')) && ...
                    (~strcmp(dl(i).name,'..'))
                run_in_subdirs([d sep dl(i).name], Suffixes, apply);
            end
        else
            [pathstr, name, ext] = fileparts(dl(i).name) ;
            for j=1: length(Suffixes)
                if ((strcmp(ext,Suffixes{j})))
                    apply([d sep dl(i).name]);
                end
            end
        end
    end
end
