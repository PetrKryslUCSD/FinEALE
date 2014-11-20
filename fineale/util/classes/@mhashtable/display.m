% displays mHashTable object contents
% *************************************************************************
% function display(mhash)
%
% Purpose: this is a mHashTable method which displays mHashTable object
%          contents
%
% Author:
%     Mikhail Poda-Razumtsev
%
% Initial Version:
%     26.11.2005
%
% Change History:
%
% *************************************************************************
function display(mhash)

isLoose = strcmp(get(0, 'FormatSpacing'), 'loose');

if(length(inputname(1)) ~= 0)
    if isLoose, disp(' '), end
    disp(sprintf('%s =', inputname(1)));
end

if isLoose, disp(' '), end

if isempty(mhash)
    fprintf('\tEmpty\n\n' );
else
    display(HashtableStore('toString', mhash.ind));
    disp(' ')
end

