function Status=pvpython (scriptfile)
    % Run the Paraview pvpython executable.
    %
    %     function Status=pvpython (scriptfile)
    %
    if ispc
        executable = 'pvpython.exe';
    elseif isunix
        executable = 'pvpython';
    else
        executable = 'pvpython';
    end
    Status=system( [executable ' ' scriptfile '']);
end
