% Format a number into a zero-padded string.
%
% function s=n2s0p(n)
% 
% Format a number into a string such that the number is zero-padded to four
% digits.
%
% n2s0p(13)
% ans =
% 0013
function s=n2s0p(n)
        if n<10
            s=['000' num2str(n)];
        elseif n<100
            s=['00' num2str(n)];
        elseif n<1000
            s=['0' num2str(n)];
        else
            s=num2str(n);
        end
    end