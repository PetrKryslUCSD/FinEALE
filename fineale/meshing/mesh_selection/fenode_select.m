function nodelist = fenode_select(fens, options)
% Select nodes.
%
% function nodelist = fenode_select(fens, options)
%
% Select nodes using some criterion, for instance and node is selected
% because it is inside a box or on its surface.
% 
% box 
% 
% Example: fenode_select(fens,struct ('box',[1 -1 2 0 4 3])), selects
% nodes which are strictly inside the box  
%  -1<= x <=1     0<= y <=2     3<= z <=4
% 
% distance
% 
% Example: fenode_select(fens,struct ('distance',0.5, 'from',[1 -1])), selects
% nodes which are Less than 0.5 units removed from the point [1 -1].
%
% nearestto
% 
% Example: fenode_select(fens,struct ('nearestto',[1 -1])), selects
% the node nearest to the point [1 -1].
%
% The option 'inflate' may be used to increase or decrease the extent of
% the box (or the distance) to make sure some nodes which would be on the
% boundary are either excluded or included.
% 
% Example: fenode_select(fens,struct ('box',[1 -1 0 0 0 0],'inflate',0.01))
% selects nodes along the line segment between x=-1, and x= 1, where all
% the nodes in the box that one gets by inflating up the segment by 0.01 in
% all directions.
% 
%
% See also: v_select

    xyz=fens.xyz;
    nodelist = v_select(xyz, options);
    nodelist =reshape(nodelist,1,length(nodelist));
end

