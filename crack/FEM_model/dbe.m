function [dbmat] = dbe(nevab,nstre,bmatx,dmatx)

format long;

for istre=1:nstre
for ievab=1:nevab
dbmat(istre,ievab)=0.0;
for jstre=1:nstre
dbmat(istre,ievab)=dbmat(istre,ievab)+dmatx(istre,jstre) ... 
	    *bmatx(jstre,ievab);
end
end
end

end %endfunction

