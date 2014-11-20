function [fens,fes,groups,edge_fes,edge_groups] =  targe2_mesher_vl(v,thickness,varargin)
    %     Automatic triangulation of the domain given by a polyline.
    %
    % function [fens,fes,groups,edge_fes,edge_groups] =
    %                          targe2_mesher_vl(v,thickness,varargin)
    %
    % Create triangulation of a domain given as an array of vertices. No holes
    % may be defined.
    %
    % Create a uniform triangulation of the domain given as an array of vertex
    % coordinates, v, one row per vertex. The vertices needs to be supplied in
    % counterclockwise order.
    %
    % Input: 
    % v= an array of vertex
    % coordinates, v, one row per vertex. The vertices needs to be supplied in
    % counterclockwise order.
    % thickness = thickness of the domain
    % varargin-- optional argument: structure with name-value pairs
    % In addition to the attributes recognized by targe2_mesher(), one can
    % specify the uniform mesh size
    %                - mesh_size:  mesh size.
    %
    % Output: The same as for targe2_mesher().
    %
    % Examples: 
    %     h=0.25;
    %     [fens,fes,groups,edge_fes,edge_groups]=...
    %         targe2_mesher_vl([-1,-1;2,-2;3,3;-1,1], 1,struct( 'mesh_size',h));
    %     gv= reset(graphic_viewer);
    %     drawmesh({fens,subset(fes,groups{1})},'gv',gv,'fes','facecolor','r','shrink', 0.9)
    %     view(2)
    %
    % See also: targe2_mesher
    

    nv =length(v);
    minSide=Inf;
    for i= 1:nv
        Commands{i}=sprintf('curve %d line %g %g %g %g\n',i,v(i,1),v(i,2),v(mod(i,nv)+1,1),v(mod(i,nv)+1,2));
        minSide=min([minSide,norm(v(i,:)-v(mod(i,nv)+1,:))]);
    end
    %     0.2
    s=sprintf('subregion 1  property 1 boundary ');
    for i= 1:length(v)
        s=[s,sprintf('%d ',i)];
    end
    Commands{end+1}=s;
    options=struct();
    if nargin >2
        options=varargin{1};
    end
    if isfield(options,'mesh_size')
        mesh_size =options.mesh_size;
    else
        mesh_size =minSide/1;
    end
    Commands{end+1}=sprintf('m-ctl-point constant %g\n',mesh_size);
    
    [fens,fes,groups,edge_fes,edge_groups] =  targe2_mesher(Commands,thickness,varargin{:});
    return;
end
