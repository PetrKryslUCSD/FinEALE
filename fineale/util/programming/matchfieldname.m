% Match fields of a struct.
%
% function out = matchfieldname(S, Name)
% 
% Match the Name on input against the names of the fields of the structure
% S. The string Name could be just a partial version, with  different
% capitalization, of the name of the field.  If any of the field names
% matches, its name is returned.
% out = matched field name or empty if there is no  field name that can be
% matched.  Note that the match is inexact.
% 
% Example:
% options = struct('initialstep', 0.01,'relTol',1e-3);
% out = matchfieldname(options, 'relT')
% out =
% relTol
%
% Copyright 2009 Petr Krysl
% 
function out = matchfieldname(S, Name)
    out = [];
    f=fieldnames(S);
    for i= 1:length(f)
        lf{i} =lower (f{i});
    end 
    j=strmatch(lower(Name),lf);
    if (~isempty(j))
        out = f{j};
    end
end
