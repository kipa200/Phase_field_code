function [gstif]= fract_stiff_v1(npoin,nelem,nnode,nstre,ndime,ndofn, ...
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

for ielem=1:nelem

 %--initialize element stiffnesses & rhs:

for ievab=1:nevab
eload(ievab) = 0.0;
for jevab=1:nevab
estif(ievab,jevab) = 0.0;
end
end

for ievab=1:nevab2
eload1(ievab) = 0.0;
for jevab=1:nevab2
estif1(ievab,jevab) = 0.0;
end
end

for ievab=1:nevab2
for jevab=1:nevab3
estif2(ievab,jevab) = 0.0;
estif3(jevab,ievab) = 0.0;
end
end

for ievab=1:nevab3
eload2(ievab) = 0.0;
for jevab=1:nevab3
estif4(ievab,jevab) = 0.0;; 
end
end

%---
%--- element nodal values     
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

[cartd,djacb,gpcod]=jacob2(ielem,elcod,kgasp,shape,deriv,nnode,ndime);

dvolu=djacb*weigp(igaus)*weigp(jgaus);

if(nnode == 3)
dvolu=djacb*weigp(igaus);
end

%-- values at the integration points:

phigp =0.0;
phirgp = 0.0;

for inode=1:nnode
phigp = phigp + ephi(inode)*shape(inode);
phirgp = phirgp + ephir(inode)*shape(inode);
end
 

mtype = 1;
[dmatx] = modps(mtype,ntype,nstre,props);

%--- 

[bmatx]=bmats(cartd,shape,inode);

[dbmat] = dbe(nevab2,nstre,bmatx,dmatx);

%element stiffness:

%--- estif1:

for ievab=1:nevab2
for jevab=1:nevab2
for istre=1:nstre

estif1(ievab,jevab)=estif1(ievab,jevab)+((1.0-phigp)^2 +constk)* ...
	    bmatx(istre,ievab)*dbmat(istre,jevab)*dvolu;

end
end
end

%---estif2:

for ievab=1:nevab2
dummy(ievab) =0.0;
end

for istre =1:nstre
for inode =1:nnode
dummy(istre,inode) =stres(ielem,kgasp,istre)*shape(inode);
end
end

for ievab=1:nevab2
for jevab=1:nevab3
for istre=1:nstre

estif2(ievab,jevab) = estif2(ievab,jevab)-2.0*(1.0-phigp)*bmatx(istre,ievab)* ...
	              dummy(istre,jevab)*dvolu;

end
end
end

%--- estif4

for ievab=1:nevab3
for jevab=1:nevab3
for istre=1:ndime
estif4(ievab,jevab) = estif4(ievab,jevab) + cenerg*constl*cartd(istre,ievab)* ...
	              cartd(istre,jevab)*dvolu;
end
end
end


%--- strain energ

senerg = 0.0;
for istre=1:nstre
senerg = senerg + 0.5*stres(ielem,kgasp,istre) * stran(ielem,kgasp,istre);	    
end

for inode=1:nnode
for jnode=1:nnode
estif4(inode,jnode) = estif4(inode,jnode) +((cenerg/constl) + 2.0*senerg)* ...
	                 shape(inode)*shape(jnode)*dvolu;

end
end

%--- penalty term:

constx = 0.0;
if(phirgp < 0.0 )
constx = -phirgp;
end

for inode=1:nnode
for jnode=1:nnode
estif4(inode,jnode) = estif4(inode,jnode) +(coneta/dtime)*constx^(constn-1)* ...
	               shape(inode)*shape(jnode)*dvolu;
end
end

end %jgaus
end %igaus

%--- global stiffness matrix:

%--- assemble estif1:

for inode =1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn2
ievab =(inode-1)*ndofn2+idofn;
itotv =(lnode-1)*ndofn2+idofn;
%
for jnode =1:nnode
knode = lnods(ielem,jnode);
for jdofn =1:ndofn2
jevab =(jnode-1)*ndofn2+jdofn;
jtotv =(knode-1)*ndofn2+jdofn;

gstif = gstif +sparse(itotv,jtotv,estif1(ievab,jevab),ntotv,ntotv);

end
end
end
end

% assemble estif2 :

for inode =1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn2
ievab =(inode-1)*ndofn2+idofn;
itotv =(lnode-1)*ndofn2+idofn;
%
for jnode =1:nnode
knode = lnods(ielem,jnode);
for jdofn =1:ndofn3
jevab =(jnode-1)*ndofn3+jdofn;
jtotv =(knode-1)*ndofn3+jdofn+ntotv2;

gstif = gstif +sparse(itotv,jtotv,estif2(ievab,jevab),ntotv,ntotv);

end
end
end
end

% assembe estif 3 as transpose of estif2:

for inode =1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn3
ievab =(inode-1)*ndofn3+idofn;
itotv =(lnode-1)*ndofn3+idofn+ntotv2;
%
for jnode =1:nnode
knode = lnods(ielem,jnode);
for jdofn =1:ndofn2
jevab =(jnode-1)*ndofn2+jdofn;
jtotv =(knode-1)*ndofn2+jdofn;

gstif = gstif +sparse(itotv,jtotv,estif2(jevab,ievab),ntotv,ntotv);

end
end
end
end

% assemble estif4:

for inode =1:nnode
lnode = lnods(ielem,inode);
for idofn =1:ndofn3
ievab =(inode-1)*ndofn3+idofn;
itotv =(lnode-1)*ndofn3+idofn+ntotv2;
%
for jnode =1:nnode
knode = lnods(ielem,jnode);
for jdofn =1:ndofn3
jevab =(jnode-1)*ndofn3+jdofn;
jtotv =(knode-1)*ndofn3+jdofn+ntotv2;

gstif = gstif +sparse(itotv,jtotv,estif4(ievab,jevab),ntotv,ntotv);

end
end
end
end

end %ielem

end %endfunction





