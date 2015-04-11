function A= extract_energy_from_abaqus_dat(File,Tagline)
 A=[];
 try
  fid=fopen(File,'r');
  while true
   temp = fgetl(fid);   if ~ischar(temp), break, end
   temp=deblank(fliplr(deblank(temp(end:-1:1))));
   if (strcmpi(temp,Tagline))
    for j=1:3
     temp = fgetl(fid);
    end
    temp = fgetl(fid);
    A = sscanf(temp(40:end), '%g');
    fclose(fid);
    return;
   end
  end
 catch,fclose(fid);end
end

