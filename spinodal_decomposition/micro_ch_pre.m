function [con] =micro_ch_pre(Nx,Ny,c0,iflag)

format long;

NxNy = Nx*Ny;

noise = 0.02;

if(iflag == 1)

for i=1:Nx
for j=1:Ny
con(i,j) =c0 + noise*(0.5-rand);
end
end

else
con=zeros(NxNy,1);

for i=1:Nx
for j=1:Ny
ii =(i-1)*Nx+j;
con(ii) =c0 + noise*(0.5-rand);
end
end

end %if

end %endfunction
