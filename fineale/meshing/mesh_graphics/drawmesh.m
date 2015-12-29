function gv = drawmesh(aMesh,varargin)
% Draw a finite element mesh.
%
% function gv = drawmesh(aMesh,varargin)
%
% Call as:
%     drawmesh(Mesh),
%     drawmesh(Mesh,'fes'), which chooses to draw just the finite elements; or
%     drawmesh(Mesh,'nodes'), which chooses to draw just the nodes; or
%     drawmesh(Mesh,'fes','nodes'), draws both;
%     drawmesh(Mesh,'fes','nodes','facecolor','none'), turns off fill
%     for element faces;
%     drawmesh(Mesh,'fes','nodes','facecolor','green','shrink', 0.8),
%     chooses green as the color with which to feel the element faces, and
%     each element is shrunk with a factor of 0.8 about its barycenter
%   where
%     Mesh is either
%        string: name of the mesh file
%             or
%        cell array of two elements, such as the output of
%              a mesh file, i.e. {fens fes}
%   Examples:
%     drawmesh('mpart','fes')
%
%     [fens fes]=H8_block(10,3,4,5,3,8);
%     mesh{1}=fens;
%     mesh{2}=fes;
%     drawmesh(mesh); %  or drawmesh({fens fes});
% 
% Optional name/value pairs (varargin):
%     gv=graphic viewer; if supplied, the graphics is added to the graphic viewer gv
%     facecolor=color to be used for the mesh faces (Default 'none');
%     edgecolor=Color to be used for the mesh edges (Default 'black');
%     linewidth=line with (default 1);
%     facealpha=face transparency (default 1, which means of opaque);
%     shrink=factor by which each element is shrunk about its barycenter (Default 1);
%
    facecolor='none';
    edgecolor='black';
    linewidth=1;
    facealpha=1;
    shrink=1;
    offset =[];
    gv=[];
    label=[];
    length_units =1;
    node_list = [];
    if (nargin == 1)
        draw_nodes=1;
        draw_fes=1;
    else
        draw_nodes=0;
        draw_fes=0;
        pos = 1;
        while (pos <= nargin-1)
            arg = varargin{pos};
            switch arg
                case 'fes'
                    draw_fes=1;
                case 'nodes'
                    draw_nodes=1;
                case 'label'
                    pos=pos+1;
                    label=varargin{pos};
                case 'facecolor'
                    pos=pos+1;
                    facecolor=varargin{pos};
                case 'edgecolor'
                    pos=pos+1;
                    edgecolor=varargin{pos};
                case 'linewidth'
                    pos=pos+1;
                    linewidth=varargin{pos};
                case 'facealpha'
                    pos=pos+1;
                    facealpha=varargin{pos};
                case 'gv'
                    pos=pos+1;
                    gv=varargin{pos};
                case 'offset'
                    pos=pos+1;
                    offset=varargin{pos};
                case 'shrink'
                    pos=pos+1;
                    shrink=varargin{pos};
                case 'node_list'
                    pos=pos+1;
                    node_list=varargin{pos};
                case 'length_units'
                    pos=pos+1;
                    length_units=varargin{pos};
                otherwise
                    error('Unrecognized argument!');
            end
            pos = pos + 1;
        end
    end
    if ~(draw_fes*draw_nodes)
        draw_fes=1;
    end

    % Mesh
    if (ischar(aMesh))
        eval(['[fens,fes] = ' aMesh ';']);
    else
        fens=aMesh{1};
        fes=aMesh{2};
    end
    % Geometry
    xyz=fens.xyz;
    geom = nodal_field(struct ('name', ['geom'], 'dim', size(xyz,2), 'fens',fens));
    limits=bounding_box(xyz);

    % Clear the figure
    if isempty(gv)
        gv=graphic_viewer;
        gv=reset (gv,struct ('limits',limits));
    end
    Rfld= nodal_field(struct ('name',['Rfield'], 'dim', 9, 'data', ones(geom.nfens,1)*reshape(eye(3),1,9)));
    context = struct ('x',geom,'u',0*geom,'R',Rfld,...
        'facecolor',facecolor,...
        'edgecolor',edgecolor,...
        'linewidth',linewidth,...
        'facealpha',facealpha,...
        'shrink', shrink,...
        'length_units', length_units);
    if (draw_fes)
        % Plot the fes
        if isempty(label)
            draw(fes,gv,context);
        else
%             for i=1:length(fes)
%                 if get(fes(i),'label')==label
%                     draw(fes(i),gv,context);
%                 end
%             end
            error(' Not implemented')
        end
    end

    if (draw_nodes)
        % Plot the node numbers
        box = [];
        if (isempty(node_list))
            node_list =1:count(fens);
        end
        for i=node_list
            box = update_box(box,xyz(i,:));
        end
        if (isempty(offset))
            offset =20*norm(diff(reshape(box,length(box)/2,2)))/(count(fens)^3);
        end
        context.offset=offset;
        context.list= node_list;
        context.length_units= length_units;
        %         context.fontsize= 18;
        draw(fens, gv, context);
    end
    % gv=interact (gv,struct ([]));

    return;
end