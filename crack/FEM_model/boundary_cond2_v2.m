function[gstif,gforce,treac] =boundary_cond2_v2(npoin,nvfix,nofix,iffix,fixed,ndofn,facto, ...
						gstif,gforce,tdisp)

format long;

ndofn2=2;
ntotv=npoin*ndofn;
treac = zeros(nvfix,ndofn2);

for ivfix = 1:nvfix
lnode = nofix(ivfix);

for idofn =1:ndofn2
if(iffix(ivfix,idofn) == 1)
itotv =(lnode - 1)*ndofn2 +idofn;

treac(ivfix,idofn) =treac(ivfix,idofn)- gstif(itotv,:)*tdisp(1:ntotv);

gstif(itotv,: ) = 0.0;

gstif(itotv,itotv) = 1.0;
gforce(itotv) = fixed(ivfix,idofn)*facto-tdisp(itotv);

end
end
end

end %endfunction
