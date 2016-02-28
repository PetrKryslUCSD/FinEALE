function Status=paraview (datafile)
    % Run the Paraview paraview executable.
    %
    % function Status=paraview (datafile)
    %
    % Paraview is a visualization tool. It can be used to visualize 
    %     results exported to a VTK file.
    %
    % See also: vtk_export_mesh
    
    if ispc 
        %         executable = 'C:\Program Files (x86)\ParaView 3.12.0\bin\paraview.exe';
        executable ='C:\Program Files\ParaView 4.4.0\bin\paraview.exe'; 
    elseif isunix
        executable = 'paraview';
    else
        executable = 'paraview';
    end
    
    if ~exist('datafile','var')
        Status=exist( [executable]);
        return;
    end
    
    if ispc 
        Status=system( ['"' executable '"' ' --data="' datafile '"&']);
    else 
        Status=system( [executable ' --data="' datafile '"&']);
    end
end
