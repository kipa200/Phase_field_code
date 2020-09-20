function [stran,stres] = stress_fract_v2(asdis,nelem,npoin,nnode,ngaus,nstre,props, ...
					 ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp, ...
					 dgdx,dvolum,tdisp)
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

stres = zeros(nelem,mgaus,nstre);
stran = zeros(nelem,mgaus,nstre);
eldis = zeros(nelem,nevab);

%--Nodal displacements:

for inode =1:nnode
lnode = lnods( :,inode);
for idofn=1:ndofn2
ievab =(inode-1)*ndofn2+idofn;	    
itotv =(lnode-1)*ndofn2+idofn;
eldis( :,ievab)=tdisp(itotv);
end
end

%--- elasticity matrix:

mtype =1;
[dmatx] = modps(mtype,ntype,nstre,props);

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

%--- strain matrix:

[bmatx]=bmats2(dgdx,nelem,nnode,nstre,nevab,kgasp);

%--- calculate the strains:

for istre =1:nstre
for ievab =1:nevab

stran( :,kgasp,istre) = stran( :,kgasp,istre) +bmatx( :,istre,ievab).*eldis( :,ievab);

end
end

%calculate stresses:

for istre =1:nstre
for jstre =1:nstre

stres( :,kgasp,istre) = stres( :,kgasp,istre) + dmatx(istre,jstre).*stran( :,kgasp,jstre);

end
end

end %igaus
end %jgaus

end %endfunction

