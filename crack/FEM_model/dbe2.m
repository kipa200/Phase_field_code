function [dbmat] = dbe2(nelem,nevab,nstre,bmatx,dmatx)

%multiply bmatx with dmatx

format long;

dbmat = zeros(nelem,nstre,nevab);

for istre=1:nstre
for ievab=1:nevab
%dbmat( :,istre,ievab)=0.0;
for jstre=1:nstre
dbmat( :,istre,ievab)=dbmat( :,istre,ievab)+dmatx(istre,jstre) ... 
            *bmatx( :,jstre,ievab);
end
end
end

end %endfunction
