function [fens,fes,groups,edge_fes,edge_groups] = targe2_mesher(fs,thickness,varargin)
% Automatic triangulation of a general 2-D domain. 
%
% function [fens,fes,groups,edge_fes,edge_groups] =
%                     targe2_mesher(fs,thickness,varargin)
%
% Automatic triangulation of a general 2-D domain. The input is given as a cell array
% of commands to the mesher Targe2. These are described in the manual,
% targe2_users_guide.pdf.
%
% Input:
% fs -- cell array of strings, each string a command recognized
%       by the mesher (targe2_users_guide.pdf) 
% thickness -- thickness of the two-dimensional slab
% varargin-- optional argument: structure with name-value pairs
%                - quadratic: convert the triangulation to
%                             quadratic elements, true or false;
%                - merge_tolerance: see this option in
%                             targe2_users_guide.pdf
%                - and any other attributes recognized by the constructors
%                of the geometric cells, for instance axisymm.
%
% Output:
%  fens= the finite element nodes
%  fes=the triangles
% The following are optional outputs:
%  groups=cell array of indexes of the triangles in the individual subregions;
%       For instance groups{3} is the list of the triangles that belong to
%       subregion 3
%  edge_fes=cell array of the edge elements
%  edge_groups=cell array of indexes of the edge elements on the individual
%       edges;
%
% Examples: 
% h=5;
% [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
%     'curve 1 line 0 0 75 0',...
%     'curve 2 line 75 0 75 40',...
%     'curve 3 line 75 40 0 0',...
%     'curve 4 circle Center 60 10 radius 8',...
%     ['subregion 1  property 1 boundary '...
%     ' 1  2 3 hole -4'],...
%     ['subregion 2  property 2 boundary 4'],...
%     ['m-ctl-point constant ' num2str(h)]
%     }, 1);
% gv= reset(graphic_viewer);
% drawmesh({fens,subset(fes,groups{1})},'gv',gv,'fes','facecolor','r','shrink', 0.9)
% drawmesh({fens,subset(fes,groups{2})},'gv',gv,'fes','facecolor','y','shrink', 0.9)
% view(2)
%
% See also: 

    in='in';
    out ='out';
    quadratic= false;
    
    if nargin >2
        options=varargin{1};
    end
    options.id=0;
    options.conn = [];
    if isfield(options,'axisymm')
        if ~options.axisymm
            options.other_dimension=thickness;
        end
    else
        options.other_dimension=thickness;
    end
    merge_tolerance=eps;
    if isfield(options,'merge_tolerance')
        merge_tolerance =options.merge_tolerance;
    else
        merge_tolerance =eps;
    end

    if nargin >2
        if isfield(options,'quadratic')
            quadratic = options.quadratic;
        end
    end

    d=pwd;
    if ~isdir(fineale_work_path)
        mkdir(fineale_work_path);
    end
    cd(fineale_work_path);
    
    if ( exist (in, 'file'))
        delete(  in); 
    end
    if ( exist (out, 'file'))
        delete(  out);
    end
    

    fid =fopen ([in], 'w');
    for i=1:length(fs)
        fprintf(fid,'%s\n',fs{i});
    end
    fclose (fid);

    if strcmp(computer,'PCWIN')
%         c=['"' finesse_path filesep 'meshing' filesep 'targe2' filesep 'Targe2.exe" -m ' num2str(options.merge_tolerance) ' -i ' in ' -f  2  -o ' out];
        c=['"' fineale_path filesep 'meshing' filesep 'triangle' filesep 'targe2' filesep 'Targe2.exe" ' ' -i ' in ' -f  2  -o ' out];
    elseif strcmp(computer,'GLNX86')  
        exec=['"' fineale_path filesep 'meshing' filesep 'triangle' filesep 'targe2' filesep 'targe2_GLNX86"'];
        c=[exec ' -i ' in ' -f  2 -o ' out ];
        system(['chmod +x ' exec]);
    elseif strcmp(computer,'GLNXA64')
        exec=['"' fineale_path filesep 'meshing' filesep 'triangle' filesep 'targe2' filesep 'targe2_GLNXA64"'];
        c=[exec ' ' ' -i ' in ' -f  2 -o ' out ];
        system(['chmod +x ' exec]);
    elseif strcmp(computer,'MACI')||strcmp(computer,'MACI64')
        exec=['"' fineale_path filesep 'meshing' filesep 'triangle' filesep 'targe2' filesep 'targe2_MACI"'];
        c=[exec ' -i ' in ' -f  2 -o ' out ];
        system(['chmod +x ' exec]);
    elseif strcmp(computer,'PCWIN64')
%         c=['"' finesse_path filesep 'meshing' filesep 'targe2' filesep 'Targe2.exe" -m ' num2str(options.merge_tolerance) ' -i ' in ' -f  2  -o ' out];
        c=['"' fineale_path filesep 'meshing' filesep 'triangle' filesep 'targe2' filesep 'Targe2.exe" ' ' -i ' in ' -f  2  -o ' out];
    else
        error(['Computer platform ' computer ' is not supported: contact the author'])
    end
    system (c);

    fid =fopen ([out], 'r');
    l=fgets(fid);
    totals =sscanf(l,'%d');
    xy =fscanf(fid,'%g', [4,totals(1)]);
    p = xy';

    es =fscanf(fid,'%g', [4,totals(2)]);
    es = es';

    ts =fscanf(fid,'%g', [5,totals(3)]);
    ts = ts';

    fclose (fid);

    %     Assign groups
    for i= 1:totals(3)
        group =ts(i,5);
        if ~exist('groups')
            groups{group} = [];
        end
        if length(groups)<group
            groups{group} = [];
        end
        groups{group} = [groups{group} i];
    end
    for i= 1:totals(2)
        group =es(i,4);
        if group~=0
            if ~exist('edge_groups')
                edge_groups{group} = [];
            end
            if length(edge_groups)<group
                edge_groups{group} = [];
            end
            edge_groups{group} = [edge_groups{group} i];
        end
    end

    nnodes= size (p, 1);
    fens=fenode_set(struct('xyz',p(p(:,1),2:3)));
    
    if quadratic
        nedges=size(es, 1);
        nmidnodes =nedges;
        es= [es, (nnodes+1:nnodes+nmidnodes)'];
        xy=fens.xyz;
        xy(nnodes+1:nnodes+nmidnodes,:)=zeros(nmidnodes,2);
        edge_conns=zeros(nedges,3);
        for i=1:nedges
            midx=p(es(i,2:3),2:3);
            xy(nnodes+i,:)=(midx(1,:) +midx(2,:))/2;
            edge_conns(i,:) =[es(i,2:3),nnodes+i];
        end
        options.conn=edge_conns;
        edge_fes=fe_set_L3(options);
        fens=fenode_set(struct('xyz',xy));
        
        ncells=size (ts, 1);
        conns=zeros(ncells,6);
        for i=1:ncells
            nns=[intersect(es(ts(i,4),2:3),es(ts(i,2),2:3)),...
                intersect(es(ts(i,2),2:3),es(ts(i,3),2:3)),...
                intersect(es(ts(i,3),2:3),es(ts(i,4),2:3))];
            enns=es(ts(i,2:4),5)';
            nns = [nns, enns];
            if det([1 p(nns(1),2:3); 1 p(nns(2),2:3); 1 p(nns(3),2:3)])<0
                nns=nns([1, 3, 2, 6, 5, 4]);
            end
            conns(i,:) =nns;
        end
        options.conn=conns;
        fes=fe_set_T6(options);
    else
        ncells=size (ts, 1);
        options.id =zeros(ncells,1);
        options.conn =zeros(ncells,3);
        for i=1:ncells
            nns=unique([es(ts(i,2),2:3) es(ts(i,3),2:3) es(ts(i,4),2:3)]);
            if det([1 p(nns(1),2:3); 1 p(nns(2),2:3); 1 p(nns(3),2:3)])<0
                nns =nns([1, 3, 2]);
            end
            options.conn(i,:) =nns;
        end
        fes=fe_set_T3(options);
        
        nedges=size(es, 1);
        options.conn =zeros(nedges,2);
        %         options.other_dimension=1.0;
        for i=1:nedges
            options.conn(i,:) =es(i,2:3);
        end
        edge_fes=fe_set_L2(options);
    end
    cd(d);

    % mesh{1}=fens;
    % mesh{2}=fes;
    % drawmesh(mesh); view (2);

    return;
end

