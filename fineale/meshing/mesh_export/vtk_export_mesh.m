function vtk_export_mesh (theFile,Connectivity,Points,Cell_types,options)
% Export mesh to a VTK 1.0 file as an unstructured grid.
%
% function vtk_export_mesh (theFile,Connectivity,Points,Cell_types,options)
%
% theFile= File name as a string,
% Connectivity= connectivity array, one row per cell,
% Points= Coordinate array, one row per point,
% Cell_types= Cell type code: L2=3, T3=5, Q4=9, T4=10, H8=12
% options = structure,  with optional attributes:
% options.binary= Boolean flag, should the file be written as binary?
% Default is false (file is ASCII).
% options.scalars =Array of per-node data, same number of rows as  array Points.
% options.scalars_name = String (no spaces), name for the scalars array.
% options.vectors =Array of per-node data, same number of rows as  array Points, three columns.
% options.vectors_name = String (no spaces), name for the vectors array.

    binary= ~true;
    scalars=[]; scalars_name= ['Data'];
    vectors=[]; vectors_name= ['Data'];
    if ( exist ( 'options', 'var') )
        if (isfield(options, 'binary'))
            binary= options.binary;
        end
        if (isfield(options, 'scalars'))
            scalars= options.scalars;
        end
        if (isfield(options, 'scalars_name'))
            scalars_name= options.scalars_name;
        end
        if (isfield(options, 'vectors'))
            vectors= options.vectors;
        end
        if (isfield(options, 'vectors_name'))
            vectors_name= options.vectors_name;
        end
    end
    if (~iscell(Connectivity))
        Connectivity={Connectivity};
    end
    if (~iscell(Cell_types))
        Cell_types={Cell_types};
    end
    [pathstr, name, ext] = fileparts(theFile);
    if (~strcmp(ext,'vtk'))
        theFile = [theFile '.vtk'];
    end
    fid=fopen(theFile,'wt');
    if (fid==-1)
        warning (['Could not open ' theFile])
        return
    end
    fprintf(fid,'# vtk DataFile Version 1.0\n');
    fprintf(fid,'FinEALE Export\n');
    if (~binary)
        fprintf(fid,'ASCII\n');
    else
        fprintf(fid,'BINARY\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS %d double\n',size(Points,1));%
    if (~binary)
        for i= 1:size(Points, 1)
            fprintf(fid,'%g %g %g\n',Points(i,:));
        end
    else
        fwrite(fid,cast(Points,'double'),'double','n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    for k=1:length(Connectivity)
        f=Connectivity{k};
        ctype=zeros(size(f, 1),1)+Cell_types{k};
        fprintf(fid,'CELLS %d %d\n',size(f,1),(size(f,1)*(size(f,2)+1)));%
        if (~binary)
            for i= 1:size(f, 1)
                fprintf(fid,'%d ',size(f,2));
                for j= 1: size(f,2)
                    fprintf(fid,'%d ',f(i,j)-1);
                end
                fprintf(fid,'\n');
            end
        else
            fwrite(fid,cast([zeros(size(f,1),1)+size(f,2),f-1],'int32'),'int32','n');
        end
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'CELL_TYPES %d\n',size(f,1));%
        if (~binary)
            for i= 1:size(f, 1)
                fprintf(fid,'%d\n',ctype(i));
            end
        else
            fwrite(fid,cast(ctype,'int32'),'int32','n');
        end
    end
    
    if (~isempty(scalars))
        fprintf(fid,'point_data %d\n',length(scalars));%
        fprintf(fid,'SCALARS %s double\n',scalars_name);%
        fprintf(fid,'LOOKUP_TABLE default\n');%
        if (~binary)
            for j= 1:length(scalars)
                    fprintf(fid,'%g\n',scalars(j));
                end
        else
            fwrite(fid,cast(scalars,'double'),'double','n');
        end
    end
    
    if (~isempty(vectors))
        fprintf(fid,'point_data %d\n',size(vectors,1));%
        fprintf(fid,'VECTORS %s double\n',vectors_name);%
        if (~binary)
            for j= 1:size(vectors,1)
                    fprintf(fid,'%g %g %g\n',vectors(j,:));
                end
        else
            fwrite(fid,cast(vectors,'double'),'double','n');
        end
    end
    
    fprintf(fid,'\n');
    fid=fclose(fid);
end