% Automatic graded triangulation of a general 2-D domain. 
%
% function [fens,fes,groups,edge_fes,edge_groups] =
%       targe2_mesher_adapt(fs,curconn,x,hests,thickness,varargin)
%
% See targe2_mesher() documentation for the input and output arguments.
% 
% Additionally:
% curconn-Current triangle connectivity  array
% hests-Desired mesh sizes at the locations x
% x-Array of locations at which mesh sizes are given. The same number of 
%      rows as there are mesh sizes in the above array, two coordinates 
%      per row.
% 
% Note: Please make sure to start the command
function [fens,fes,groups,edge_fes,edge_groups] = ...
    targe2_mesher_adapt(fs,curconn,x,hests,thickness,varargin)
    s=fs{1}(1:6);
    if (~strcmp(upper(s),'REGEXT'))
        fs=cat(2,{['REGEXT ' num2str(min(x)) ' ' num2str(max(x))]},fs);
    end
    for j=1:size(curconn,1)
        conn=curconn(j,:);
        s=['mc-t ' num2str(x(conn(1),:)) '  ' num2str(x(conn(2),:)) ' ' num2str(x(conn(3),:)) ' ' num2str(min(hests(conn))) ' ' num2str(min(hests(conn))) ' ' num2str(min(hests(conn)))];
        fs{end+1}=s;
    end 
    [fens,fes,groups,edge_fes,edge_groups] = targe2_mesher(fs,thickness,varargin{:});
end

