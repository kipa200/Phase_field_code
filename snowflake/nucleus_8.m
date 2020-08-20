function [phi,tempr] = nucleus(Nx,Ny,seed)

format long;

%for i=1:Nx
%for j=1:Ny

%phi(i,j) = 0.0;
%tempr(i,j) = 0.0;

%end
%end 

phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);

for i=1:Nx
for j=1:Ny
if ((i-Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-2*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-3*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-4*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-5*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-6*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-7*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
if ((i-8*Nx/8)^2+(j-Ny/2)^2 < seed)
  phi(i,j) = 1.0;
end
end
end

end %endfunction
