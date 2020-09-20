function [stran,stres] = stress_fract_v1(asdis,nelem,npoin,nnode,ngaus,nstre,props, ...
				 ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp, ...
				 tdisp)
format long;

%--- order of integration
ngaus2=ngaus;
if(nnode == 3)
ngaus2=1;
end

mgaus = ngaus*ngaus2;

%--initialize:

ndofn2 =2;
nevab2 =nnode*ndofn;
ntotv2 =npoin*ndofn2;
nevab = nnode*ndofn2;

for ielem =1:nelem

for igaus =1:mgaus
for istre =1:nstre
stres(ielem,igaus,istre) = 0.0;
stran(ielem,igaus,istre) = 0.0;
end
end

%--Nodal displacements:

for inode =1:nnode
lnode = lnods(ielem,inode);
for idofn=1:ndofn2
ievab =(inode-1)*ndofn2+idofn;	    
itotv =(lnode-1)*ndofn2+idofn;
eldis(ievab)=tdisp(itotv);
end
end

%--- elasticity matrix:

mtype =1;
[dmatx] = modps(mtype,ntype,nstre,props);

%--- Coordinates of element nodes:

for inode=1:nnode
lnode=lnods(ielem,inode);
for idime=1:ndime
elcod(idime,inode)=coord(lnode,idime);
end
end

%--- calculate strains and stresses at integration points:

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

%--- strain matrix:

[bmatx]=bmats(cartd,shape,inode);

%--- calculate the strains:

for istre =1:nstre
for ievab =1:nevab

stran(ielem,kgasp,istre) = stran(ielem,kgasp,istre) +bmatx(istre,ievab)*eldis(ievab);

end
end

%calculate stresses:

for istre =1:nstre
for jstre =1:nstre

stres(ielem,kgasp,istre) = stres(ielem,kgasp,istre) + dmatx(istre,jstre)*stran(ielem,kgasp,jstre);

end
end

end %igaus
end %jgaus

end % ielem
end %endfunction

