function [fens,fes]=Comsol_mesh_import(filename,xyzscale,saveReadData)
% Import mesh from the Comsol mesh file.
%
% function [fens,fes]=Comsol_mesh_import(filename,xyzscale,saveReadData) 
%
% Import tetrahedral (4- and 10-node) Comsol Mesh (.mphtxt).
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
        
        % find the number of nodes
        temp = '';
        while 1
            temp = fgetl(fid); 
            if (temp==-1)
                error('No ''number of mesh points'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'numberofmeshpoints'))
                nnode=sscanf(temp, '%d');
                break;
            end
        end
        % Allocate the node array
        node=zeros(nnode,3);
        
        % find the lowest mesh point index
        lowest=0;
        temp = '';
        while 1
            temp = fgetl(fid); 
            if (temp==-1)
                error('No ''lowest mesh point index'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'lowestmeshpointindex'))
                lowest=sscanf(temp, '%d');
                break;
            end
        end
        
        % Find the mesh point coordinates
        temp = '';
        while 1
            temp = fgetl(fid);   
            if (temp==-1)
                error('No ''mesh point coordinates'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'Meshpointcoordinates'))
                break;
            end
        end
        
        % read the nodes
        for j=1:nnode
            temp = fgetl(fid); 
            if (temp==-1)
                error('Ran out of data for nodes')
            end
            A = sscanf(temp, '%g %g %g');
            node(j,:)=A;
        end 
        
        % Find the solid element block
        temp = '';
        while 1
            temp = fgetl(fid);   
            if (temp==-1)
                error('No ''type name'' heading for 3 tet in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'typename'))
                if (regexpi(temp,'\s*3\s*tet\s*#.*'))
                    break;
                end
            end
        end
        
        % find the number of nodes per element
        ennod=0;
        temp = '';
        while 1
            temp = fgetl(fid);
            if (temp==-1)
                error('No ''number of nodes per element'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'numberofnodesperelement'))
                ennod=sscanf(temp, '%d');
                break;
            end
        end
        
        % find the number of elements
        nelem=0;
        temp = '';
        while 1
            temp = fgetl(fid);
            if (temp==-1)
                error('No ''number of elements'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'numberofelements'))
                nelem=sscanf(temp, '%d');
                break;
            end
        end
        % Allocate the element array
        element=zeros(nelem,ennod);
        
        % Find the elements
        temp = '';
        while 1
            temp = fgetl(fid);   
            if (temp==-1)
                error('No ''elements'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'elements'))
                break;
            end
        end
        
        % Now read elements
        for j=1:nelem
            temp = fgetl(fid) ;
            if (temp==-1)
                error('Ran out of data for elements')
            end
            A = sscanf(temp, '%d', ennod);
            element(j,:)=A';
        end 
        
        % Adjust the index into the node array
        if (lowest <1)% In Matlab the indexes are one-based
            element=element+1;
        end
        
        % find the number of geometric entity indices
        gein=0;
        temp = '';
        while 1
            temp = fgetl(fid) ;
            if (temp==-1)
                error('No ''number of geometric entity indices'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'numberofgeometricentityindices'))
                gein=sscanf(temp, '%d');
                break;
            end
        end
        % Allocate the element array
        gei=zeros(gein,1);
        
        
        % Find the geometric entity indices
        temp = '';
        while 1
            temp = fgetl(fid);   
            if (temp==-1)
                error('No ''geometric entity indices'' heading in file')
            end
            [n]=heading_name(temp);
            if (strcmpi(n,'Geometricentityindices'))
                break;
            end
        end
        
        % Now read geometric entity indices
        for j=1:nelem
            temp = fgetl(fid); 
            if (temp==-1)
                error('Ran out of data for indices')
            end
            A = sscanf(temp, '%d', 1);
            gei(j,:)=A';
        end 
        
        
        fclose(fid);
    end
    
    % Create output arguments. First the nodes
    clear fens
    xyz =xyzscale*node;
    fens=fenode_set(struct('xyz',xyz));
    
    % Now the geometric cells for each element set
    clear fes
    switch ennod
        case 4
            Constructor =@feset_T4;
        case 10
            Constructor =@feset_T10;
        otherwise
            error('Unknown element type');
    end
    % Number of different element sets
    elsets= unique(gei);
    for j=1:length(elsets)
        ix= find(gei==elsets(j));
        fes{j}= Constructor(struct('conn',element(ix,:),'label',elsets(j)));
    end
    
    % If there was just one element set, return a single feset object
    if (length(elsets)==1)
        fes=fes{1};
    end
    
    % Should the data that we read be saved?
    if (saveReadData && (Doread))
        save(saveFile,'xyz','elsets')
    end
end
 
function [n]=heading_name(str)
n=[];
f=strfind(str, '#');
if (isempty(f))
    return;
end
n=str(f+1:end);
c=1;
for j=1:length(n)
    if (n(j)~=' ')
        n(c)=n(j);
        c=c+1;
    end
end
n=n(1:c-1);
end
