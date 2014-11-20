function clean_sweep(Suffixes)
    % Clean up by removing  designated files.
%
% function clean_sweep(Suffixes)
%
% Clean up in the current directory and its subdirectories by removing designated files.
%
% Suffixes = cell array of suffixes of files to remove (default is 
%           {'.asv','.log','.bak','.m~1~','.m~2~','.m~3~'})
%
% See also: sniff_out

    if  ~exist('Suffixes')
        Suffixes ={'.asv','.log','.bak','.m~1~','.m~2~','.m~3~'};
    end
    function   remove_found(x)
        %         x
        try
            delete([x]);
        catch
            disp([x ' could not be deleted'])
        end
    end
%     sniff_out(Suffixes)
    sniff_out(Suffixes,@remove_found)
    return;
end
