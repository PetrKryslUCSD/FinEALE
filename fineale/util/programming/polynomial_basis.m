function N= polynomial_basis(ndimensions, order, p)
% Calculate point value of polynomial basis.
%
% function N= polynomial_basis(ndimensions, order, p)
%
% Example: Symbolic evaluation.
% polynomial_basis(3, 1, [sym('x'),sym('y'),sym('z')])

N=zeros((order+1)^ndimensions,1,class(p));
switch ndimensions
    case 3
        t=1;
        for i=0:1:order
            for j=0:1:order
                for k=0:1:order
                    N(t)=p(1)^i*p(2)^j*p(3)^k; t=t+1;
                end
            end
        end
    otherwise
        N=[];
end
N=reshape(N,[],1);


 
         