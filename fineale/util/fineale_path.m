% Return the path to the fineale topmost folder.
%
% function p = fineale_path
%
function p = fineale_path
    s='fineale_init';
    p = which(s);
    p = p(1:length(p)-length(s)-3);
end
