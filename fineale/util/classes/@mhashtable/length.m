% object isempty overloading
% *************************************************************************
% function flag = length(obj)
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
function flag = length(obj)

flag = HashtableStore('size', obj.ind);
