% Extract just a few elements from an array.
% 
% function out = somel(A,varargin)
% 
% The function returns elements of the array given by the indexes supplied
% as varargin.
%
% Example:
% somel(abs(diag(D)),[2,3]) %: give me just the second and third number
function out = somel(A,varargin)
    out=A(varargin{:});
end