function [posgp, weigp] =gauss(ngaus,nnode)

format long;

if(nnode == 3)

if(ngaus == 1)
posgp(1)=1.0/3.0;
posgp(2)=1.0/3.0;
weigp(1)=0.5;
end

if(ngaus == 3)
posgp(1)=0.5;
posgp(2)=0.5;
posgp(3)=0.0;
%
posgp(4)=0.0;
posgp(5)=0.5;
posgp(6)=0.5;

weigp(1)=1.0/6.0;
weigp(2)=1.0/6.0;
weigp(3)=1.0/6.0;
end

if(ngaus == 7)
posgp(1)=0.0;
posgp(2)=0.5;
posgp(3)=1.0;
posgp(4)=0.5;
posgp(5)=0.0;
posgp(6)=0.0;
posgp(7)=1.0/3.0;

posgp(8)=0.0;
posgp(9)=0.0;
posgp(10)=0.0;
posgp(11)=0.5;
posgp(12)=1.0;
posgp(13)=0.5;
posgp(14)=1.0/3.0;

weigp(1)=1.0/40.0;
weigp(2)=1.0/15.0;
weigp(3)=1.0/40.0;
weigp(4)=1.0/15.0;
weigp(5)=1.0/40.0;
weigp(6)=1.0/15.0;
weigp(7)=9.0/40.0
end

end



if (nnode ~=3)

if(ngaus == 2)
posgp(1)=-0.57735026918963;
weigp(1)=1.0;
end
if(ngaus > 2)
posgp(1)=-0.7745966241483;
posgp(2)=0.0;
weigp(1)=0.55555555555556;
weigp(2)=0.88888888888889;
end

kgaus=ngaus/2;
for igash=1:kgaus
jgash=ngaus+1-igash;
posgp(jgash)=-posgp(igash);
weigp(jgash)=weigp(igash);
end
end

end %endfunction
