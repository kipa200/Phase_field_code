function[tdisp,stres,stran] =initilize(nelem,npoin,nnode,ngaus,ndofn,...
						     nstre,nvfix,isolve,icase)
  
format long;

ndofn2 = ndofn-1;
ntotv2 = npoin*ndofn2;
ntotv = npoin*ndofn;

ngaus2 = ngaus;
if(nnode == 3)
ngaus2=1;
end
mgaus = ngaus*ngaus2;

%-- initialize stress and strain:

if(isolve == 1)

for ielem =1:nelem
for igaus=1:mgaus
for istre=1:nstre
 
stres(ielem,igaus,istre) = 0.0;
stran(ielem,igaus,istre) = 0.0;
end
end
end

for itotv=1:ntotv
tdisp(itotv,1)=0.0;
end

end %if

if(isolve == 2)

stres = zeros(nelem,mgaus,nstre);
stran = zeros(nelem,mgaus,nstre);

tdisp =zeros(ntotv,1);

end %if

if(icase == 1)

ncrack =30;
crack(1) = 2010;
crack(2) = 2011;
crack(3) = 2012;
crack(4) = 2013;
crack(5) = 2014;
crack(6) = 2015;
crack(7) = 2016;
crack(8) = 2017;
crack(9) = 2018;
crack(10) = 2019;
crack(11) = 2020;
crack(12) = 2021;
crack(13) = 2022;
crack(14) = 2023;
crack(15) = 2024;
%
crack(16) = 2051;
crack(17) = 2052;
crack(18) = 2053;
crack(19) = 2054;
crack(20) = 2055;
crack(21) = 2056;
crack(22) = 2057;
crack(23) = 2058;
crack(24) = 2059;
crack(25) = 2060;
crack(26) = 2061;
crack(27) = 2062;
crack(28) = 2063;
crack(29) = 2064;
crack(30) = 2065;

for icrack=1:ncrack
lnode = crack(icrack);
itotv=ntotv2+lnode;
tdisp(itotv) = 0.99;
end

end

if(icase == 2) 

ncrack =26;
crack(1) = 291;
crack(2) = 292;
crack(3) = 296;
crack(4) = 297;
crack(5) = 1553;
crack(6) = 1554;
crack(7) = 1555;
crack(8) = 1556;
crack(9) = 1557;
crack(10) = 1558;
crack(11) = 1559;
crack(12) = 1560;
crack(13) = 1561;
crack(14) = 1562;
crack(15) = 1563;
crack(16) = 1568;
crack(17) = 1569;
crack(18) = 1570;
crack(19) = 1571;
crack(20) = 1572;
crack(21) = 1573;
crack(22) = 1574;
crack(23) = 1575;
crack(24) = 1576;
crack(25) = 1577;
crack(26) = 1578;

for icrack=1:ncrack
lnode = crack(icrack);
itotv=ntotv2+lnode;
tdisp(itotv) = 0.99;
end

end 


end %endfunction

