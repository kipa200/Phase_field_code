%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                    %       
%          FEM PHASE-FIELD CODE FOR FRACTURE         %
%                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get initial wall time:
time0=clock();

icase =1;

global in;
if(icase == 1)
in = fopen('fract_1ca.inp','r');
end

if(icase == 2)
in = fopen('fract_1hca.inp','r');
end

global out;
out = fopen('result_1.out','w');

global out2;
out2 = fopen('force-disp','w');

%--- input time integration parameters:

nstep=    4000; 
nprnt=      25;
dtime=1.0;

miter= 5;
toler= 5.0e-3;

tfacto = 0.0;
dfacto = 0.0005;

isolve =2;

% Material spefici parameters:

constk = 1.0e-6;
cenerg = 0.001;
constl = 0.125;   
constn = 2;
coneta = 20.0e5*2;

%-----------
%input data:
%-----------

[npoin,nelem,nvfix,ntype,nnode,ndofn,ndime,ngaus, ...
      nstre,nmats,nprop,lnods,matno,coord,props,nofix,  ...
      iffix,fixed]=input_fem_elast();

ntotv = npoin*ndofn;
ndofn2 = 2;
ntotv2 = npoin*ndofn2;

[tdisp,stres,stran] =initilize(nelem,npoin,nnode,ngaus,ndofn,...
					     nstre,nvfix,isolve,icase);

[posgp, weigp] =gauss(ngaus,nnode);

if(isolve == 2)

[dgdx,dvolum] = cart_deriv(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			   ngaus,ntype,lnods,coord,posgp,weigp);
end

%--- time integration

for istep = 1:nstep

tfacto = tfacto + dfacto;

tdisp_old = tdisp;
gforce = zeros(ntotv,1);

if(isolve == 1)

[gstif]= fract_stiff_v1(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			ngaus,ntype,lnods,coord,props,posgp, ...
			weigp,dtime,constk,cenerg,constl,constn, ...
			coneta,stres,stran,tdisp,tdisp_old);

end
%
if(isolve == 2)

[gstif]= fract_stiff_v2(npoin,nelem,nnode,nstre,ndime,ndofn, ...
			ngaus,ntype,lnods,coord,props,posgp, ...
			weigp,dgdx,dvolum,dtime,constk,cenerg, ...
			constl,constn,coneta,stres,stran,tdisp, ...
			tdisp_old);

end

%--- newton iteration:

for iter = 1:miter

%--- boundary conditions:

if(isolve == 1)

[gstif,gforce,treac] =boundary_cond2_v1(npoin,nvfix,nofix,iffix,fixed,ndofn, ...
					tfacto,gstif,gforce,tdisp);
end

if(isolve == 2)

[gstif,gforce,treac] =boundary_cond2_v2(npoin,nvfix,nofix,iffix,fixed,ndofn,...
					tfacto,gstif,gforce,tdisp);
end

%--- solve:

asdis = gstif\gforce;

%--- update:

tdisp = tdisp + asdis;

%--- adjust small deviations:

if(isolve == 1)
for ipoin=1:npoin
itotv = ntotv2 +ipoin;
if(tdisp(itotv) > 0.999)
  tdisp(itotv) = 0.999;
end
if(tdisp(itotv) < 0.0)
  tdisp(itotv) = 0.0;
end
end

end 
%
if(isolve == 2)
dummy =tdisp(ntotv2+1:ntotv);

inrange = (dummy > 0.999);
dummy(inrange) =1.0;
inrange = ( dummy < 0.0);
dummy(inrange) = 0.0;
tdisp(ntotv2+1:ntotv) = dummy;
end

%---calculate stress & strain increments:

if (isolve == 1)

[stran,stres] = stress_fract_v1(asdis,nelem,npoin,nnode,ngaus,nstre,props, ...
				ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp, ...
				tdisp);
end

if(isolve == 2)
[stran,stres] = stress_fract_v2(asdis,nelem,npoin,nnode,ngaus,nstre,props, ...
				ntype,ndofn,ndime,lnods,matno,coord,posgp,weigp, ...
				dgdx,dvolum,tdisp);
end

%--- check norm for convergence:

normF = norm(gforce,2);

if(normF <= toler)
break;
end

%--- calculate residual force vector:

if(isolve == 1)

[gforce]= residual_v1(npoin,nelem,nnode,nstre,ndime,ndofn, ...
		      ngaus,ntype,lnods,coord,props,posgp, ...
		      weigp,dtime,constk,cenerg,constl,constn, ...
		      coneta,stres,stran,tdisp,tdisp_old);
end

if(isolve ==2 )

[gforce]= residual_v2(npoin,nelem,nnode,nstre,ndime,ndofn, ...
		      ngaus,ntype,lnods,coord,props,posgp, ...
		      weigp,dgdx,dvolum,dtime,constk,cenerg, ...
		      constl,constn,coneta,stres,stran,tdisp, ...
		      tdisp_old);
end

end %end of Newton

%--- print data for force-disp curves:

lnode = nofix(nvfix);

nvfix2 =nvfix/2;

sumr = 0.0;
for ivfix=1:nvfix2
sumr = sumr +treac(ivfix,2);
end

fprintf(out2,'%14.6e %14.6e\n',tdisp((lnode-1)*2+2),sumr);

%--- print results:

if(mod(istep,nprnt) == 0 )

fprintf('Done step: %5d\n',istep);

%fname=sprintf('time_%d.out',istep);
%out=fopen(fname,'w');
%ntotv2 =npoin*ndofn2;
%for ipoin = 1:npoin
%itotv =ntotv2 + ipoin;
%fprintf(out,'%14.6e %14.6e %14.6e\n',coord(ipoin,1), coord(ipoin,2),tdisp(itotv));
%end
%fclose(out);

%--- write to vtk file with updated mesh

for ipoin=1:npoin
for idofn=1:ndofn2
itotv =(ipoin-1)*ndofn2+idofn;	    
cord2(ipoin,idofn) =coord(ipoin,idofn) + 10.0*tdisp(itotv);
jtotv=ntotv2+ipoin;
cont1(ipoin)=tdisp(jtotv);
end
end

write_vtk_fem(npoin,nelem,nnode,lnods,cord2,istep,cont1);

end %if

end %istep

compute_time = etime(clock(),time0)

fprintf(out,'compute time: %7d\n',compute_time);

