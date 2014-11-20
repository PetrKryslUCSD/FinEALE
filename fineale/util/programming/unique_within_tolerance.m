% Find unique values in array within tolerance.
%
% function xo=unique_within_tolerance(x,xtolerance)
%
% unique_within_tolerance works as unique, but with a numerical tolerance
% for deciding which numbersin the array are unique. 
% 
% C = unique_within_tolerance(A,eps) for the array A, returns the same
% values as in A except for those that are within machine epsilon of each
% other.   The values  that are within the tolerance are merged using the
% mean as the output value. The values of C are in sorted order.
%
% See also: unique

function xo=unique_within_tolerance(x,xtolerance)
    x= unique(x);
    xd  =diff(x);
    xo=x;
    i= find(xd<=xtolerance);
    if ( ~isempty(i))
         for ii=1:length(i)
             mx=mean([x(i(ii)),x(i(ii)+1)]);
            x(i(ii))=mx;,x(i(ii)+1)=mx;
         end
         xo=x(setdiff(1:length(x),i));
    end
end
       