function [dgdx,dvolum] = cart_deriv(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			    ngaus,ntype,lnods,coord,posgp,weigp)


ngaus2 = ngaus;
if(nnode == 3)
ngaus2 = 1;
end

mgaus = ngaus*ngaus2;

dvolum = zeros(nelem,mgaus);

dgdx = zeros(nelem,mgaus,ndime,nnode);

for ielem = 1:nelem			     
for inode=1:nnode
lnode=lnods(ielem,inode);
for idime=1:ndime
elcod(idime,inode)=coord(lnode,idime);
end
end

%=== gauss points:

kgasp=0;
for igaus=1:ngaus
exisp=posgp(igaus);
for jgaus=1:ngaus2
etasp =posgp(jgaus);
if(nnode ==3)
etasp=posgp(ngaus+igaus);
end


kgasp=kgasp+1;
mgaus=mgaus+1;

[shape,deriv]=sfr2(exisp,etasp,nnode);
[cartd,djacb,gpcod]=jacob2(ielem,elcod,kgasp,shape,deriv,nnode,ndime);

dvolu=djacb*weigp(igaus)*weigp(jgaus);

if(nnode == 3)
dvolu=djacb*weigp(igaus);
end

dvolum(ielem,kgasp)=dvolu;

for idime=1:ndime
for inode=1:nnode
dgdx(ielem,kgasp,idime,inode)=cartd(idime,inode);
end
end

end %igaus
end %jgaus
end % ielem

end

%endfunction




