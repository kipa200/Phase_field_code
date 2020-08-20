function [grad] =laplacian(nx,ny,dx,dy)

format long;

nxny=nx*ny;

r=zeros(1,nx);
r(1:2)=[2,-1];
T=toeplitz(r);

E=speye(nx);

grad=-(kron(T,E)+kron(E,T));

%-- for periodic boundaries

for i=1:nx
ii=(i-1)*nx+1;
jj=ii+nx-1;
grad(ii,jj)=1.0;
grad(jj,ii)=1.0;

kk=nxny-nx+i;
grad(i,kk)=1.0;
grad(kk,i)=1.0;
end

grad = grad /(dx*dy);

end %endfunction

