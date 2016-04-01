function Status=abaqus_job(inpfile)
    % Run the Abaqus executable.
    %
    % function Status=abaqus_job(inpfile)
    %
    % The Abaqus executable is proprietary to SIMULIA. If you have access
    % to it, it is possible to run it from within FinEALE solution files.
    % One may have to input the complete path to the executable. See below.
    
    if ispc 
        executable = 'C:\SIMULIA\Abaqus\6.14-2SE\code\bin\abq6142se.exe';
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
