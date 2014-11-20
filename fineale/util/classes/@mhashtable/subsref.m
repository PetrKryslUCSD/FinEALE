% manages Java-like methods syntax
% *************************************************************************
% function varargout = subsref(obj, index)
%
% Purpose: the methods syntax for mHashTable class is not Matlab-like but
%          Java-like. This is achieved using subsref
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
function ArgOut = subsref(obj, index)

% check calling mode
if index(1).type ~= '.', return; end

% get action and arguments
arg     = {};
action  = index(1).subs;
if length(index) > 1, arg = index(2).subs; end

% always return one argument
ArgOut = HashtableStore(action, obj.ind, arg{:});
