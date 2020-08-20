function [phi,tempr] = nucleus(Nx,Ny,seed)

format long;

phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);

for i=1:Nx
for j=1:Ny
if ((i-Nx/2)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
end
end

end %endfunction
