%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%    PHASE-FIELD FINITE-DIFFIRENCE CODE FOR SOLVING           %
%                  CAHN-HILLIARD EQUATION                     %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%== get intial wall time:
time0=clock();
format long;

out2 = fopen('time_energ.out','w')

%-- Simulation cell parameters:

Nx = 64;
Ny = 64;
NxNy= Nx*Ny;
dx = 1.0;
dy = 1.0; 

%--- Time integration parameters:

nstep =   10000; %时间步总数
nprint=      50; %输出频率
dtime =  1.0e-2; %时间积分增量     
ttime =     0.0; %总时间

%--- Material specific Parameters:

c0 = 0.40; %流动系数
mobility = 1.0; %梯度能量系数
grad_coef= 0.5;

%
%---prepare microstructure:
%

iflag =1; %Generate initial grain microstructure.iflag = 1 for bi-crystal and iflag = 2 is forpolycrystal.
[con] = micro_ch_pre(Nx,Ny,c0,iflag);

%
%--- Evolve
%

for istep =1:nstep

ttime = ttime +dtime;

for i=1:Nx
for j=1:Ny

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

if(im == 0)
 im=Nx;
end
if(ip == (Nx+1))
  ip=1;
end

if(jm == 0) 
  jm = Ny;
end

if(jp == (Ny+1))
  jp=1;
end

hne=con(ip,j);
hnw=con(im,j);
hns=con(i,jm);
hnn=con(i,jp);
hnc=con(i,j);

lap_con(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);

%--- derivative of free energy:

[dfdcon]=free_energ_ch_v1(i,j,con);

dummy(i,j) = dfdcon - grad_coef*lap_con(i,j);

end %for j
end %for i

%--

for i=1:Nx
for j=1:Ny

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

if(im == 0)
  im=Nx;
end
if(ip == (Nx+1))
  ip=1;
end

if(jm == 0) 
  jm = Ny;
end

if(jp == (Ny+1))
  jp=1;
end

hne=dummy(ip,j);
hnw=dummy(im,j);
hns=dummy(i,jm);
hnn=dummy(i,jp);
hnc=dummy(i,j);

lap_dummy(i,j) = (hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);

%-- time integration:

con(i,j) = con(i,j) + dtime*mobility*lap_dummy(i,j);

%-- for small deviations:

if(con(i,j) >= 0.9999);
con(i,j)= 0.9999;
end
if(con(i,j) < 0.00001);
con(i,j) = 0.00001;
end

end %j
end %i

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

%-- calculate total energy

[energ] = calculate_energ(Nx,Ny,con,grad_coef);

fprintf(out2,'%14.6e %14.6e\n',ttime,energ);

write_vtk_grid_values(Nx,Ny,dx,dy,istep,con);  

end %if 

end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);

