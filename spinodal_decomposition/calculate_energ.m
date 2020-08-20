function [energ] = calculate_energ(Nx,Ny,con,grad_coef)

format long;

energ =0.0;

for i=1:Nx-1
ip = i + 1;
for j=1:Ny-1 
jp = j + 1;

energ = energ + con(i,j)^2*(1.0-con(i,j))^2 + ...
  0.5*grad_coef*((con(ip,j)-con(i,j))^2 + (con(i,jp)-con(i,j))^2);
end
end

end %endfunction

