% object isempty overloading
% *************************************************************************
% function flag = isempty(obj)
%
% Purpose: 
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
function flag = isempty(obj)

flag = HashtableStore('isempty', obj.ind);
