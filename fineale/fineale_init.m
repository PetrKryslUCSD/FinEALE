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
    % set the homepath to the directory containing thei file
    homePath = fileparts(which(mfilename()));
    addpath(fullfile (homePath, 'util'));
    
    if ( ~exist('options','var')) || (isempty(options ))
        options.useCXSparse= true;
    end
    
    % Set paths
    addpath(homePath,...
        fullfile (homePath, 'classes'),...
        fullfile (homePath, 'meshing'), ...
        fullfile (homePath, 'algorithms'),...
        fullfile (homePath, 'examples'),...
        fullfile (homePath, 'documents'),...
        fullfile (homePath, 'tutorials'),...
        fullfile (homePath, 'util'));
    add_subdirs_to_path (fullfile (homePath, 'classes'), sep);
    add_subdirs_to_path (fullfile (homePath, 'meshing'), sep);
    add_subdirs_to_path (fullfile (homePath, 'util'), sep);
    add_subdirs_to_path (fullfile (homePath, 'algorithms'), sep);
    add_subdirs_to_path (fullfile (homePath, 'examples'), sep);
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
    	if (exist(fullfile(homePath, '..', 'CXSparse'), 'dir'))
    	    addpath(fullfile(homePath, '..', 'CXSparse', 'MATLAB', 'CSparse'),....
    	    	    fullfile(homePath, '..', 'CXSparse', 'MATLAB', 'UFget'));
    	end
    end
    
    %     if (usejava('desktop') )
    %         doc fineale;
    %     end
    
    return;

end


function add_subdirs_to_path(d,sep)
    dl=dir(d);
    for i=1:length(dl)
        if (dl(i).isdir)
            if      (~strcmp(dl(i).name,'.')) && ...
                    (~strcmp(dl(i).name,'..')) && ...
                    (~strcmp(dl(i).name,'CVS')) && ...
                    (~strcmp(dl(i).name,'cvs')) && ...
                    (~strcmp(dl(i).name,'private')) && ...
                    (~strcmp(dl(i).name(1),'@'))
                addpath(fullfile(d,dl(i).name))
                add_subdirs_to_path(fullfile(d,dl(i).name), sep);
            end
        end
    end
    return;
end
 