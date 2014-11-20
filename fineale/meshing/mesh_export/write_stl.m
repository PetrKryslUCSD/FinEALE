function write_stl (theFile,f,v)
    % Write out  surface mesh in the STL format.
    %
    % function write_stl (theFile,f,v)
    %
    % Arguments
    % theFile= filename,
    % f= connectivity array, one row per polygon,
    % v= vertex location array, one per vertex

if (size(f,2)~=3)
    if (size(f,2)~=4)
        error (['Expected triangles or Quadrilaterals']);
    else
        f=[f(:,[1,2,3]);f(:,[3,4,1])];
    end
end
    [pathstr, name, ext] = fileparts(theFile);
    if (~strcmp(ext,'stl'))
        theFile = [theFile '.stl'];
    end
    fid=fopen(theFile,'w');
    if (fid==-1)
        warning (['Could not open ' theFile])
        return
    end
    fprintf(fid,'solid %s\n',name);
    for i= 1:size(f, 1)
        a=v(f(i,2),:)-v(f(i,1),:);
        b=v(f(i,3),:)-v(f(i,1),:);
        fprintf(fid,'facet normal %g %g %g\n',skewmat(a)*b');
        fprintf(fid,'outer loop\n');
        for j= 1:3
            fprintf(fid,'vertex %g %g %g\n',v(f(i,j),:));
        end
        fprintf(fid,'endloop\n');
        fprintf(fid,'endfacet\n');
    end
    fprintf(fid,'endsolid %s\n',name);
    fid=fclose(fid);
end