function Status=abaqus_job(inpfile)
    % Run the Abaqus executable.
    %
    % function Status=abaqus_job(inpfile)
    %
    %
    
    if ispc 
        executable = 'abaqus';
    elseif isunix
        executable = 'abaqus';
    else
        executable = 'abaqus';
    end
    
    if ~exist('inpfile','var')
        Status=exist( [executable]);
        return;
    end
    
    if (~exist( [executable]))
        error(['Executable "'  executable '" not found'])
    end
    
    if ispc 
        Status=system( ['"' executable '"' ' job="' inpfile '"']);
    else 
        Status=system( [executable ' job="' inpfile '"']);
    end
end
