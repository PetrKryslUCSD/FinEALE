% object isempty overloading
% *************************************************************************
% function [m, n] = size(obj)
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
function [m, n] = size(obj)

m = HashtableStore('size', obj.ind);
if m == 0, n = 0; else n = 1; end 
