% Simple illustration of permutations.
clear p a b ip 
p=randperm(5)
a=rand(1,5)
b=a(p)
ip(p)=(1:5)
b(ip)