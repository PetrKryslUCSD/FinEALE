function [varargout] = adeal(x)
% Deal values from input array to a number of outputs
% 
% function [varargout] = adeal(x)
% 
%
% Examples: 
% [a,b,c,d]  =adeal([1,2,3,4])
    nout = max(nargout,1);
    s = length(x);
    for k=1:nout, varargout{k}= x(k); end
end
