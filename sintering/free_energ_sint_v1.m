function [dfdcon,dfdeta] =free_energ_sint_v1(i,j,con,eta,etas,npart,iflag)

format long;

A=16.0;
B= 1.0;

dfdcon =0.0;
dfdeta =0.0;


if(iflag == 1)
%--
%-- 濃度場自由能迭代:
%--

sum2 = 0.0;
sum3 = 0.0;

for ipart = 1:npart

sum2 = sum2 + etas(i,j,ipart).^2;
sum3 = sum3 + etas(i,j,ipart).^3;
end

dfdcon = B*(2.0*con(i,j) + 4.0*sum3 - 6.0*sum2) - 2.0*A*con(i,j)^2 .* ...
        (1.0-con(i,j)) + 2.0*A*con(i,j)*(1.0-con(i,j))^2;

end %if 

if(iflag == 2)
%--
%-- 相場自由能迭代:
%--

sum2 =0.0;

for ipart = 1: npart

sum2 = sum2 + etas(i,j,ipart)^2;
end

dfdeta = B*(-12.0*eta(i,j)^2 .* (2.0-con(i,j)) +12.0 *eta(i,j)* (1.0-con(i,j)) + ...
	    12.0*eta(i,j)*sum2);
end %if

end  %endfunction
