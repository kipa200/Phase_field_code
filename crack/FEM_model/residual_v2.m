function [rload ]= residual_v2(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			       ngaus,ntype,lnods,coord,props,posgp, ...
			       weigp,dgdx,dvolum,dtime,constk,cenerg, ...
			       constl,constn,coneta,stres,stran,tdisp, ...
			       tdisp_old)

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

%--- element loads:

eload =  zeros(nelem,nevab);
eload1 = zeros(nelem,nevab2);
eload2 = zeros(nelem,nevab3);

%--
%--- elemental values     
%--

ephi = zeros(nelem,nnode);
ephir = zeros(nelem,nnode);

for inode =1:nnode
lnode = lnods( :,inode);
itotv = ntotv2 +lnode;
ephi( :,inode) = tdisp(itotv);
ephir( :,inode) = tdisp(itotv)-tdisp_old(itotv);
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

%--- values at the  integration points:

phigp =zeros(nelem,1);
phirgp =zeros(nelem,1);

for inode=1:nnode
phigp = phigp + ephi( :,inode)*shape(inode);
phirgp = phirgp + ephir( :,inode)*shape(inode);
end

mtype = 1;
[dmatx] = modps(mtype,ntype,nstre,props);

%--- cartesien derivative matrices;

[bmatx1]=bmats1(dgdx,nelem,nnode,nstre,nevab2,kgasp);
[bmatx2]=bmats2(dgdx,nelem,nnode,nstre,nevab2,kgasp);


% strain energy:

senerg =zeros(nelem,1);
for istre=1:nstre

senerg = senerg + 0.5*stres( :,kgasp,istre) .* stran( :,kgasp,istre);	    

end

%--- penalty term:

constx = zeros(nelem,1);
inrange = (phirgp < 0 );
constx(inrange) = -phirgp(inrange);

%--- residuals (rhs)

%--- eload1:

for ievab=1:nevab2
for istre=1:nstre

eload1( :,ievab) = eload1( :,ievab) + ((1.0-phigp).^2 + constk) .*bmatx2( :,istre,ievab).* ...
	    stres( :,kgasp,istre).*dvolum( :,kgasp);

end
end

%--- eload2

dummy = zeros(nelem,ndime);

for istre = 1:ndime
for ievab = 1:nevab3
dummy( :,istre) = dummy( :,istre) + bmatx1( :,istre,ievab).*phigp;
end
end

for ievab=1:nevab3
for istre=1:ndime
eload2( :,ievab) = eload2( :,ievab) + cenerg*constl*bmatx1( :,istre,ievab).* ...
	    dummy( :,istre).*dvolum( :,kgasp);
end
end

for inode=1:nnode
eload2( :,inode) = eload2( :,inode) + ((cenerg/constl)+2.0*senerg).* ...
	    phigp*shape(inode).*dvolum( :,kgasp);
end

for inode=1:nnode
eload2( :,inode) = eload2( :,inode) - 2.0*shape(inode)* ... 
	    (senerg-0.5*(coneta/dtime)*constx.^constn).*dvolum( :,kgasp);

end

end %jgaus
end %igaus

%--- assemble global residuals:

for inode=1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn2
ievab=(inode-1)*ndofn2+idofn;
itotv=(lnode-1)*ndofn2+idofn;
rload(itotv) = rload(itotv) +eload1( :,ievab);
end
end

for inode=1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn3
ievab=(inode-1)*ndofn3+idofn;
itotv=(lnode-1)*ndofn3+idofn+ntotv2;
rload(itotv) = rload(itotv) + eload2( :,ievab);
end
end

end %endfunction





