utips  = [];
for j=length(us).*[15, 30, 45, 60]/60
    u=us{j};
            utip=(gather_values(u,cncl));
    utips= [   utips,     -utip(3)];
end
utips