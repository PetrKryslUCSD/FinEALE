% is this object a valid isMHashTable
% *************************************************************************
% function flag = isMHashTable(obj)
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
function flag = isMHashTable(obj)

flag = HashtableStore('isMHashTable', obj.ind);
