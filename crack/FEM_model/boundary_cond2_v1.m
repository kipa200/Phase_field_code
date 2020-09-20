function[gstif,gforce,treac] =boundary_cond2_v1(npoin,nvfix,nofix,iffix,fixed,ndofn,facto, ...
						gstif,gforce,tdisp)


global out1;

format long;

ndofn2=2;
ntotv=npoin*ndofn;

for ivfix=1:nvfix
for idofn=1:ndofn2
treac(ivfix,idofn) =0.0;
end
end

for ivfix = 1:nvfix
lnode = nofix(ivfix);

for idofn =1:ndofn2
if(iffix(ivfix,idofn) == 1)
itotv =(lnode - 1)*ndofn2 +idofn;
for jtotv = 1:ntotv

treac(ivfix,idofn) =treac(ivfix,idofn)- gstif(itotv,jtotv)*tdisp(jtotv);
gstif(itotv,jtotv) = 0.0;

end

gstif(itotv,itotv) = 1.0;
gforce(itotv) = fixed(ivfix,idofn)*facto-tdisp(itotv);
end
end
end

end %endfunction
