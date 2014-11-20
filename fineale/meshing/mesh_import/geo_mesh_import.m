function [fens,fes]=geo_mesh_import(filename,xyzscale,saveReadData)
% Import mesh from the Cosmos mesh file.
%
% function [fens,fes]= geo_mesh_import(filename,xyzscale,saveReadData) 
%
% Import tetrahedral (4- and 10-node) Cosmos Mesh (.geo).
% Limitations: only the *NODE and *ELEMENT  sections are read. Only 4-node and 10-node
% tetrahedra  are handled.
%
% Arguments:
% filename= filename,
% xyzscale= scale the coordinates by this number when they are read in,
% saveReadData= should the data of the saved in a format that allows 
%           for speedy reading of the same data later on?
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% See also: mesh_import
%    
    
    if (~exist('xyzscale','var'))
        xyzscale = 1.0;
    end
    
    if (~exist('saveReadData','var'))
        saveReadData = ~true;
    end
    
    [pathstr, name, ext] = fileparts(filename) ;
    if (isempty(pathstr)),pathstr ='.';end
    saveFile = [pathstr filesep name ,'_xyzconn' '.mat'];
    
    Doread = true;
    if (saveReadData)
        if (exist(saveFile,'file'))
            try,
                load(saveFile);
                ennod =size(conn,2);
                Doread= false;
            catch,end;
        end
    end
    if (Doread)
        fid=fopen(filename,'r');
        
        if fid == -1,
            error('Unable to open specified file.');
        end
        
        chunk=5000;
        
        nnode=0;
        node=zeros(chunk,4);
        nelem=0;
        elem=zeros(chunk,12);
        ennod=[];
        temp = '';
        while true
            temp = fgetl(fid);
            if ~ischar(temp),
                More_data=~true;
                break,
            end
            temp=strtrim(temp);
            if (strncmpi(temp,'ND',2))
                nnode=nnode+1;
                temp=strrep (temp,',',' ');
                A = sscanf(temp, 'ND %g %g %g %g');% ND, 135,      0.34087241      0.15976152    -0.038224775
                node(nnode,:)=A;
                if (size(node,1)<=nnode)
                    node(size(node,1)+1:size(node,1)+chunk,:)=zeros(chunk,4);
                end
            elseif (strncmpi(temp,'EL',2))
                nelem=nelem+1;
                temp=strrep (temp,',',' ');
                A=Parse_EL_line(temp);
                elem(nelem,1:length(A))=A;
                if (size(elem,1)<=nelem)
                    elem(size(elem,1)+1:size(elem,1)+chunk,:)=zeros(chunk,12);
                end
            end
        end
        
        node=node(1:nnode,:);
        elem=elem(1:nelem,:);
        
        fclose(fid);
        
        if (norm((1:nnode)'-node(1:nnode,1)))
            error('Nodes are not in serial order');
        end
        
        % Process output arguments
        % Extract coordinates
        xyz=node(node(1:nnode,1),2:4);
        xyz =xyzscale*xyz;
        % Cleanup element connectivities
        ennod=unique(elem(:,2));
        if  length(ennod)~=1
            error('This function cannot treat a mixture of element types at this point');
        end
        conn=elem(:,3:3+ennod-1);
    end
    
    % Create output arguments. First the nodes
    clear fens
    fens=fenode_set(struct('xyz',xyz));
    
    % Now the geometric cells for each element set
    clear fes
    switch ennod
        case 4
            Constructor =@fe_set_T4;
        case 10
            Constructor =@fe_set_T10;
        otherwise
        error('Unknown element type');
    end
    fes= Constructor(struct('conn',conn(:,1:ennod),'label',j));
    
    % Should the data that we read be saved?
    if (saveReadData && (Doread))
        save(saveFile,'xyz','elsets')
    end
end

function A=Parse_EL_line(String)
    A=[];
    if ~ischar(String), return; end
    tokens= all_tokens(String, ' ');
    if (strcmpi(tokens{1},'EL'))
        A(1)=sscanf(tokens{2}, '%g');% Element number
        A(2)=sscanf(tokens{5}, '%g');% nodes per element, 4 or 10
        for k=1:A(2)
            A(2+k)=sscanf(tokens{5+k}, '%g');% Node number
        end 
    end
end
function tokens= all_tokens(String, Separator)
    tokens={};
    while ~isempty( String )
        [tokens{end+1}, String]=strtok(String,Separator);
    end
end
