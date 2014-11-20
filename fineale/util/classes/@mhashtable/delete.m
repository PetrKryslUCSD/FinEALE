% deletes the hashtable object form memory, manual garbage collection
% *************************************************************************
% function varargout = delete(obj)
%
% Purpose: frees the memory for the object which is not used any more
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
function success = delete(obj)

% always return one argument
success = HashtableStore('delete', obj.ind);
