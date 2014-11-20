function run_all_examples(dirl)
% Run all in examples in the current directory and its subdirectories.
% 
% function run_all_examples(dirl)
% 
% dirl== optional list of directories (folders) as a cell array
%
% The results are written into a log file.
% Check that at the bottom of the log file where the script reports success
% or failure; warnings will be presented for those examples for
% which the run failed.
% 
    if (~ exist('dirl' , 'var' ))
        dirl={pwd};
    end
    if strcmp(version('-release'),'13')
        fid = fopen(['results-' '.log'],'w');
    else
        fid = fopen(['results' '-on-' datestr(now, 'mm-dd-yyyy') '-at-' datestr(now, 'HH-MM') '.log'],'w');
    end
    failed={};
    for i= 1:length(dirl)
        failed=run_in_subdirs(dirl{i}, fid, failed);
    end 
    
    fprintf(fid,'\n\n\n%s\n',['Summary:']);
    if (isempty(failed))
        fprintf(fid,'%s\n',['Completed successfully']);
    else
        for i=1:length(failed )
            fprintf(fid,'%s\n',['Warning: Failed in ' failed{i} '!!!']);
        end 
    end
    fclose (fid);
    return;
end

function failed=run_in_subdirs(d, fid, failed)
    sep=filesep;
    dl=dir(d);
    for i=1:length(dl)
        if (dl(i).isdir)
            if      (~strcmp(dl(i).name,'.')) & ...
                    (~strcmp(dl(i).name,'..')) & ...
                    (~strcmp(dl(i).name,'CVS')) & ...
                    (~strcmp(dl(i).name,'cvs')) & ...
                    (~strcmp(dl(i).name(1),'@'))
                failed=run_in_subdirs([d sep dl(i).name], fid, failed);
            end
        else
            if ((~strcmp(dl(i).name,'Contents.m')) &...
                    (~strcmp(dl(i).name,'run_all_examples.m')) &...
                    strcmp(dl(i).name,...
                    [dl(i).name(1:length(dl(i).name)-2) '.m']))
                fn=dl(i).name(1:length(dl(i).name)-2);
                disp(['Running ' fn]);
                wd=pwd;
                cd(d);
                fprintf(fid,'%s',['Running ' fn ':']);
                clear fineale_test_passed
                try
                    eval(fn);
                    if (exist ('fineale_test_passed'))
                        if (fineale_test_passed)
                            fprintf(fid,' %s.\n',['completed, results verified']);
                            disp(['         Done, verified: ' fn]);
                        else
                            failed{end+1}=fn;
                            fprintf(fid,' %s.\n',['failed, results not verified']);
                            disp(['         Failed, not verified: ' fn]);
                        end
                    else
                        fprintf(fid,' %s.\n',['completed']);
                        disp(['         Done: ' fn]);
                    end
                catch
                    failed{end+1}=fn;
                    fprintf(fid,' %s.\n',['failed to complete']);
                    disp(['         Failed: ' fn]);
                end
                cd(wd);
                pause(1)
                while ~isempty(get(0,'CurrentFigure'))
                    try
                        delete(get(0,'CurrentFigure'));
                    catch
                    end
                end
            end
        end
    end
    return;
end