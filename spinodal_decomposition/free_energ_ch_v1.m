function [dfdcon] =free_energ_ch_v1(i,j,con)

format long;

A=1.0;

dfdcon =A*(2.0*con(i,j)*(1-con(i,j))^2-2.0*con(i,j)^2 * (1.0-con(i,j)));

end %endfunction
