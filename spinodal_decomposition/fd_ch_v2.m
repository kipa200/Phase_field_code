%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%    PHASE-FIELD FINITE-DIFFIRENCE CODE FOR SOLVING           %
%                  CAHN-HILLIARD EQUATION                     %
%             (OPTIMIZED FOR MATLAB/OCTAVE)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%== get intial wall time:
time0=clock();
format long;

out2 = fopen('time_energ.out','w');

%-- Simulation cell parameters:

Nx = 64;
Ny = 64;
NxNy= Nx*Ny;
dx = 1.0;
dy = 1.0; 

%--- Time integration parameters:

nstep =    20000;
nprint=      100;
dtime =   1.0e-2;     
ttime =      0.0;

%--- Material specific Parameters:

c0 = 0.40;
mobility = 1.0;
grad_coef= 0.5;

%
%---prepare microstructure:
%

iflag =2;

[con] = micro_ch_pre(Nx,Ny,c0,iflag);

%--- Get Laplacian templet:

[grad] =laplacian(Nx,Ny,dx,dy);

%
%--- Evolve:
%

for istep =1:nstep

ttime = ttime+dtime;

%--- derivative of free energy:

[dfdcon]=free_energ_ch_v2(Nx,Ny,con);

lap_con = grad*con;

lap_con2 =grad*(dfdcon - grad_coef*lap_con);

%--- Time integration:

con = con + dtime * mobility * lap_con2;

%-- for small deviations:

inrange =(con >= 0.9999);
con(inrange) = 0.9999;
inrange =(con < 0.00001);
con(inrange) = 0.00001;

%---- print results

if((mod(istep,nprint) == 0) || (istep == 1) )

fprintf('done step: %5d\n',istep);

%fname1 =sprintf('time_%d.out',istep);
%out1 = fopen(fname1,'w');

%for i=1:Nx
%for j=1:Ny
%fprintf(out1,'%5d %5d %14.6e\n',i,j,con(i,j));
%end
%end

%fclose(out1);

%--- write vtk file:

for i=1:Nx
for j=1:Ny
ii=(i-1)*Nx+j;
con2(i,j) = con(ii);
end
end

%-- calculate total energy

[energ] = calculate_energ(Nx,Ny,con2,grad_coef);

fprintf(out2,'%14.6e %14.6e\n',ttime,energ);

write_vtk_grid_values(Nx,Ny,dx,dy,istep,con2);  

end %if 

end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);

