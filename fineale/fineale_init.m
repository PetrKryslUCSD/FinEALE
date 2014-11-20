% Initialize the fineale toolbox.
%
% Usage: fineale_init (i.e.  without arguments) or fineale_init(options)
% Note: Run in the top fineale folder.
%
% This functions purpose is to initialize the toolbox before its first use
% in a Matlab session.  Matlab needs to have appropriate subfolders of 
% fineale to its path selected can find classes and auxiliary functions.
%
% options==  optional argument, struct with  optional fields
%     useCXSparse=  should we use CXSparse?   True or false.
%
% Example: fineale_init(struct('useCXSparse',1))
%
function fineale_init(options)
    sep = filesep;
    addpath(['.' sep 'util']);
    homePath = fineale_path;
    
    if ( ~exist('options','var')) || (isempty(options ))
        options.useCXSparse= true;
    end
    
    % Set paths
    addpath(homePath,...
        [homePath sep 'classes'],...
        [homePath sep 'meshing'], ...
        [homePath sep 'algorithms'],...
        [homePath sep 'examples'],...
        [homePath sep 'documents'],...
        [homePath sep 'tutorials'],...
        [homePath sep 'util']);
    add_subdirs_to_path ([homePath sep 'classes'], sep);
    add_subdirs_to_path ([homePath sep 'meshing'], sep);
    add_subdirs_to_path ([homePath sep 'util'], sep);
    add_subdirs_to_path ([homePath sep 'algorithms'], sep);
    add_subdirs_to_path ([homePath sep 'examples'], sep);
    if strcmp(version('-release'),'13')
        warning off MATLAB:m_warning_end_without_block
    end

    disp(' ');
    disp(['fineale (C) 2014, Petr Krysl.']);
    disp(['   See the README  file for details.']);
    disp('   -------------------------------------------------------------');
    disp('   To get started, run some of the scripts in the "examples" ')
    disp('   directory. Help is available: type doc fineale (or help fineale).')
    disp(' ');
    
    Have_CXSparse = ~ (strcmp('',which ('cs_lsolve')));
    if (Have_CXSparse)
        try
            cs_lsolve(sparse([1]),[1]);
        catch
            Have_CXSparse = false ;
        end
    end
    if (~Have_CXSparse) && exist('options','var') && isfield( options,'useCXSparse') && ( options.useCXSparse )
    	if (exist([homePath sep '..' sep 'CXSparse']))
    	    addpath([homePath sep '..' sep 'CXSparse/MATLAB/CSparse'],....
    	    	    [homePath sep '..' sep 'CXSparse/MATLAB/UFget']);
    	end
    end
    
    %     if (usejava('desktop') )
    %         doc fineale;
    %     end
    
    return;

    function add_subdirs_to_path(d, sep)
        dl=dir(d);
        for i=1:length(dl)
            if (dl(i).isdir)
                if      (~strcmp(dl(i).name,'.')) & ...
                        (~strcmp(dl(i).name,'..')) & ...
                        (~strcmp(dl(i).name,'CVS')) & ...
                        (~strcmp(dl(i).name,'cvs')) & ...
                        (~strcmp(dl(i).name,'private')) & ...
                        (~strcmp(dl(i).name(1),'@'))
                    addpath([d sep dl(i).name])
                    add_subdirs_to_path([d sep dl(i).name], sep);
                end
            end
        end
        return;
    end
end
 