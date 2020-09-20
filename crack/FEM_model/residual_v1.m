function [rload ]= residual_v1(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			       ngaus,ntype,lnods,coord,props,posgp, ...
			       weigp,dtime,constk,cenerg,constl,constn, ...
			       coneta,stres,stran,tdisp,tdisp_old)

format long;

%--- order of integration

ngaus2=ngaus;
if(nnode == 3)
ngaus2=1;
end

mgaus = ngaus*ngaus2;

%--initialize global and local rhs vectors:

ntotv = npoin*ndofn;
nevab = nnode*ndofn;
ndofn2 = ndofn-1;
ndofn3 = ndofn-2;
ntotv2 = npoin*ndofn2;
nevab2 = nnode*ndofn2;
nevab3 = nnode*ndofn3;

%---global residuals:

rload = zeros(ntotv,1);

for ielem=1:nelem

%--- initialize element loads:

for ievab=1:nevab
eload(ievab) = 0.0;
end

for ievab=1:nevab2
eload1(ievab) = 0.0;
end

for ievab=1:nevab3
eload2(ievab) =0.0;
end

%--
%--- elemental values     
%--

for inode =1:nnode
lnode = lnods(ielem,inode);
itotv = ntotv2 +lnode;
ephi(inode) = tdisp(itotv);
ephir(inode) = tdisp(itotv)-tdisp_old(itotv);
end

%--- Coordinates of element nodes:

for inode=1:nnode
lnode=lnods(ielem,inode);
for idime=1:ndime
elcod(idime,inode)=coord(lnode,idime);
end
end

%--- integrate elemental loads:  

kgasp=0;
for igaus=1:ngaus
exisp=posgp(igaus);
for jgaus=1:ngaus2
etasp =posgp(jgaus);
if(nnode ==3)
etasp=posgp(ngaus+igaus);
end

kgasp=kgasp+1;

[shape,deriv]=sfr2(exisp,etasp,nnode);

[cartd,djacb,gpcod]=jacob2(ielem,elcod,kgasp,shape,deriv,nnode,ndime);

dvolu=djacb*weigp(igaus)*weigp(jgaus);

if(nnode == 3)
dvolu=djacb*weigp(igaus);
end

%--- values at the  integration points:

phigp = 0.0;
phirgp = 0.0;

for inode=1:nnode
phigp = phigp + ephi(inode)*shape(inode);
phirgp = phirgp + ephir(inode)*shape(inode);
end

mtype = 1;
[dmatx] = modps(mtype,ntype,nstre,props);

%---;

[bmatx]=bmats(cartd,shape,inode);

% strain energy:

senerg = 0.0;
for istre=1:nstre

senerg = senerg + 0.5*stres(ielem,kgasp,istre) * stran(ielem,kgasp,istre);	    

end

%--- penalty term:

constx = 0.0;
if(phirgp < 0.0 )
constx = -phirgp;
end

%--- residuals (rhs)

%--- eload1:

for ievab=1:nevab2
for istre=1:nstre

eload1(ievab) = eload1(ievab) + ((1.0-phigp)^2 + constk)*bmatx(istre,ievab)* ...
	        stres(ielem,kgasp,istre)*dvolu;

end
end

%--- eload2

for istre=1:ndime
dummy(istre) = 0.0;
end

for istre = 1:ndime
for ievab = 1:nevab3
dummy(istre) = dummy(istre) + cartd(istre,ievab)*phigp;
end
end

for ievab=1:nevab3
for istre=1:ndime
eload2(ievab) = eload2(ievab) + cenerg*constl*cartd(istre,ievab)* ...
	        dummy(istre)*dvolu;
end
end

for inode=1:nnode
eload2(inode) = eload2(inode) + ((cenerg/constl)+2.0*senerg)* ...
	        phigp*shape(inode)*dvolu;
end

for inode=1:nnode
eload2(inode) = eload2(inode) - 2.0*shape(inode)* ... 
	        (senerg-0.5*(coneta/dtime)*constx^constn)*dvolu;

end

end %jgaus
end %igaus

%--- assemble global residuals:

for inode=1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn2
ievab=(inode-1)*ndofn2+idofn;
itotv=(lnode-1)*ndofn2+idofn;
rload(itotv) = rload(itotv) +eload1( :,ievab);
end
end

for inode=1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn3
ievab=(inode-1)*ndofn3+idofn;
itotv=(lnode-1)*ndofn3+idofn+ntotv2;
rload(itotv) = rload(itotv) + eload2( :,ievab);
end
end

end %ielem

end %endfunction





