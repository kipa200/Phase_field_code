function [dfdcon] =free_energ_ch_v2(Nx,Ny,con)

format long;

A=1.0;

dfdcon =A*(2.0*con .* (1-con).^2 -2.0*con.^2 .* (1.0-con));

end %endfunction
