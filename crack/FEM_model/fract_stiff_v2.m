function [gstif]= fract_stiff_v2(npoin,nelem,nnode,nstre,ndime,ndofn, ...
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

%--initialize global and local stiffness:

ntotv = npoin*ndofn;
nevab = nnode*ndofn;
ndofn2 = ndofn-1;
ndofn3 = ndofn-2;
ntotv2 = npoin*ndofn2;
nevab2 = nnode*ndofn2;
nevab3 = nnode*ndofn3;

%---global stiffness:

gstif = sparse(ntotv,ntotv);

%--- element stiffnesses

estif = zeros(nelem,nevab,nevab);
estif1 = zeros(nelem,nevab2,nevab2);
estif2 = zeros(nelem,nevab2,nevab3);
estif3 = zeros(nelem,nevab3,nevab2);
estif4 = zeros(nelem,nevab3,nevab3); 

eload = zeros(nelem,nevab);
eload1 = zeros(nelem,nevab2);
eload2 = zeros(nelem,nevab3);

%---
%--- element nodal values     
%--

ephi = zeros(nelem,nnode);
ephir = zeros(nelem,nnode);

for inode =1:nnode
lnode = lnods( :,inode);
itotv = ntotv2 +lnode;
ephi( :,inode) = tdisp(itotv);
ephir( :,inode) = tdisp(itotv)-tdisp_old(itotv);
end

%--- integrate element stiffness  

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

%-- values at the integration points:

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

[dbmat] = dbe2(nelem,nevab2,nstre,bmatx2,dmatx);

%element stiffness:

%--- estif1:

for ievab=1:nevab2
for jevab=1:nevab2
for istre=1:nstre

estif1( :,ievab,jevab)=estif1( :,ievab,jevab)+((1.0-phigp).^2 +constk).* ...
	    bmatx2( :,istre,ievab).*dbmat( :,istre,jevab).*dvolum(:,kgasp);

end
end
end

%---estif2:

dummy = zeros(nelem,nevab2);

for istre =1:nstre
for inode =1:nnode
dummy( :,istre,inode) =stres( :,kgasp,istre)*shape(inode);
end
end

for ievab=1:nevab2
for jevab=1:nevab3
for istre=1:nstre

estif2( :,ievab,jevab) = estif2( :,ievab,jevab)-2.0*(1.0-phigp).*bmatx2( :,istre,ievab).* ...
	    dummy( :,istre,jevab).*dvolum( :,kgasp);

end
end
end

%--- estif4

for ievab=1:nevab3
for jevab=1:nevab3
for istre=1:ndime
estif4( :,ievab,jevab) = estif4( :,ievab,jevab) + cenerg*constl*bmatx1( :,istre,ievab).* ...
	    bmatx1( :,istre,jevab).*dvolum( :,kgasp);
end
end
end

%--- strain energ

senerg =zeros(nelem,1);
for istre=1:nstre
senerg = senerg + 0.5*stres( :,kgasp,istre) .* stran( :,kgasp,istre);	    
end

for inode=1:nnode
for jnode=1:nnode
estif4( :,inode,jnode) = estif4( :,inode,jnode) +((cenerg/constl) + 2.0*senerg).* ...
	    shape(inode)*shape(jnode).*dvolum( :,kgasp);

end
end


%--- penalty term:

constx = zeros(nelem,1);
inrange = (phirgp < 0 );
constx(inrange) = -phirgp(inrange);


for inode=1:nnode
for jnode=1:nnode
estif4( :,inode,jnode) = estif4( :,inode,jnode) +(coneta/dtime)*constx.^(constn-1)* ...
	                 shape(inode)*shape(jnode).*dvolum( :,kgasp);
end
end

end %jgaus
end %igaus

%--- global stiffness matrix

%--- assemble estif1:

for inode =1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn2
ievab =(inode-1)*ndofn2+idofn;
itotv =(lnode-1)*ndofn2+idofn;
%
for jnode =1:nnode
knode = lnods( :,jnode);
for jdofn =1:ndofn2
jevab =(jnode-1)*ndofn2+jdofn;
jtotv =(knode-1)*ndofn2+jdofn;

gstif = gstif +sparse(itotv,jtotv,estif1( :,ievab,jevab),ntotv,ntotv);

end
end
end
end

% assemble estif2 :

for inode =1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn2
ievab =(inode-1)*ndofn2+idofn;
itotv =(lnode-1)*ndofn2+idofn;
%
for jnode =1:nnode
knode = lnods( :,jnode);
for jdofn =1:ndofn3
jevab =(jnode-1)*ndofn3+jdofn;
jtotv =(knode-1)*ndofn3+jdofn+ntotv2;

gstif = gstif +sparse(itotv,jtotv,estif2( :,ievab,jevab),ntotv,ntotv);

end
end
end
end

% assembe estif 3 as transpose of estif2:

for inode =1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn3
ievab =(inode-1)*ndofn3+idofn;
itotv =(lnode-1)*ndofn3+idofn+ntotv2;
%
for jnode =1:nnode
knode = lnods( :,jnode);
for jdofn =1:ndofn2
jevab =(jnode-1)*ndofn2+jdofn;
jtotv =(knode-1)*ndofn2+jdofn;

gstif = gstif +sparse(itotv,jtotv,estif2( :,jevab,ievab),ntotv,ntotv);

end
end
end
end

% assemble estif4:

for inode =1:nnode
lnode = lnods( :,inode);
for idofn =1:ndofn3
ievab =(inode-1)*ndofn3+idofn;
itotv =(lnode-1)*ndofn3+idofn+ntotv2;
%
for jnode =1:nnode
knode = lnods( :,jnode);
for jdofn =1:ndofn3
jevab =(jnode-1)*ndofn3+jdofn;
jtotv =(knode-1)*ndofn3+jdofn+ntotv2;

gstif = gstif +sparse(itotv,jtotv,estif4( :,ievab,jevab),ntotv,ntotv);

end
end
end
end

end %endfunction





