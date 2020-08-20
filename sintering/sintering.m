%== 获得初始时间，设置数据显示格式:
time0=clock();
format long;

%-- 设置单元尺寸:

Nx = 100;
Ny = 100;
NxNy= Nx*Ny;

dx = 0.5;
dy = 0.5;

%--- 设置积分步:

nstep =      5000;
nprint =       50;
dtime=  1.0e-4;

%--- 材料相关属性:
npart = 9; %粒子数2或9
%
coefm = 5.0;
coefk = 2.0;
coefl = 5.0;
%
dvol = 0.040;
dvap = 0.002;
dsur = 16.0;
dgrb = 1.6;
%--准备结果存储空间:

X = 0:0.5:49.5;
Y = 0:0.5:49.5;
Zcon = zeros(100, 100);
Zeta1 = zeros(100, 100);
Zeta2 = zeros(100, 100);
%----
%微结构初始化,程序设置了粒子数为2，9两种情况:
%----
for ipart =1:npart
    for i=1:Nx
        for j=1:Ny
            con(i,j) =0.0;
            etas(i,j,ipart) = 0.0;
        end
    end
end % ipart
%---设置9粒子圆心
if(npart ~= 2)
    
    R = 10.0;
    
    xc(1)=29.0;
    yc(1)=50.0;
    
    xc(2)=50.0;
    yc(2)=50.0;
    
    xc(3)=71.0;
    yc(3)=50.0;
    
    xc(4)=50.0;
    yc(4)=29.0;
    
    xc(5)=50.0;
    yc(5)=71.0;
    
    xc(6)=39.0;
    yc(6)=39.0;
    
    xc(7)=61.0;
    yc(7)=39.0;
    
    xc(8)=39.0;
    yc(8)=61.0;
    
    xc(9)=61.0;
    yc(9)=61.0;
    
    for ipart=1:npart
        
        Rx = R;
        
        if(ipart > 5 )
            Rx = 0.5*R;
        end
        
        for i=1:Nx
            for j=1:Ny
                
                xx1 =sqrt((i-xc(ipart))^2 +(j-yc(ipart))^2);
                
                if( xx1 <= Rx)
                    con(i,j)= 0.999;
                    etas(i,j,ipart) =0.999;
                    
                end %if
            end % j
        end % i
    end % ipart
end % if

%---
if(npart == 2)
    
    R1 = 20.0;
    R2 = 0.5*R1;
    
    x1 = Nx/2;
    y1 = 40.0;
    y2 = 70.0;
    
    for i=1:Nx
        for j=1:Ny
            
            xx1 =sqrt((i-x1)^2 +(j-y1)^2);
            xx2 =sqrt((i-x1)^2 +(j-y2)^2);
            
            if( xx1 <= R1)
                con(i,j)=0.999;
                etas(i,j,1) =0.999;
            end % if
            
            if(xx2 <= R2)
                con(i,j) = 0.999;
                etas(i,j,1)=0.0;
                etas(i,j,2)=0.9999;
            end % if
        end % j
    end % i
end %if

%---
%-- 初始化 eta:
for i=1:Nx
    for j=1:Ny
        eta(i,j) = 0.0;
    end
end

for istep =1:nstep
    
    %-- 浓度场迭代:
    iflag = 1;
    for i=1:Nx
        for j=1:Ny
            
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
            
            %--- 计算自由能对浓度的导数:
            
            [dfdcon,dfdeta]=free_energ_sint_v1(i,j,con,eta,etas,npart,iflag);
            dummy(i,j) = dfdcon - 0.5*coefm*lap_con(i,j);
            
        end %for j
    end %for i
    
    %--
    
    for i=1:Nx
        for j=1:Ny
            
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
            
            lap_dummy(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
            
            %-- Mobility插值函数:
            
            phi = con(i,j)^3 *(10.0-15.0*con(i,j) + 6.0*con(i,j)^2);
            
            sum =0.0;
            for ipart =1:npart
                for jpart =1:npart
                    if(ipart ~= jpart)
                        sum = sum + etas(i,j,ipart)*etas(i,j,jpart);
                    end
                end
            end
            
            mobil = dvol*phi + dvap*(1.0-phi) + dsur*con(i,j)*(1.0-con(i,j)) + dgrb*sum;
            
            %-- 时间积分:
            
            con(i,j) = con(i,j) + dtime*mobil*lap_dummy(i,j);
            
            %-- 处理偏差值:
            
            if(con(i,j) >= 0.9999);
                con(i,j)= 0.9999;
            end
            
            if(con(i,j) < 0.00001);
                con(i,j) = 0.00001;
            end
            Zcon(i,j) = con(i,j);
        end %j
    end %i
    
    %-- 相场迭代:
    iflag = 2;
    for ipart = 1:npart
        
        for i=1:Nx
            for j=1:Ny
                eta(i,j) = etas(i,j,ipart);
            end
        end
        
        for i=1:Nx
            for j=1:Ny
                
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
                
                hne=eta(ip,j);
                hnw=eta(im,j);
                hns=eta(i,jm);
                hnn=eta(i,jp);
                hnc=eta(i,j);
                
                lap_eta(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
                
                %--- 自由能对相的导数:
                [dfdcon,dfdeta]=free_energ_sint_v1(i,j,con,eta,etas,npart,iflag);
                %-- 时间积分
                
                eta(i,j) = eta(i,j) - dtime * coefl*(dfdeta - 0.5 *coefk*lap_eta(i,j));
                
                %-- 偏差值处理:
                
                if(eta(i,j) >= 0.9999);
                    eta(i,j) = 0.9999;
                end
                
                if(eta(i,j) < 0.0001);
                    eta(i,j) = 0.0001;
                end
            end %j
        end %i
        
        %--
        
        for i=1:Nx
            for j=1:Ny
                etas(i,j,ipart) = eta (i,j);
                if(ipart==1)
                    Zeta1(i, j) =  etas(i,j,ipart);
                else
                    Zeta2(i, j) =  etas(i,j,ipart);
                end
            end%end for
        end
    end %ipart
    
    %---- 结果输出
    
    if((mod(istep,nprint) == 0) || (istep == 1) )
        subplot(1,3,1)
        pcolor(X, Y, Zcon)
        axis([0 50 0 50]);
        
        subplot(1,3,2)
        pcolor(X, Y, Zeta1)
        axis([0 50 0 50]);
        
        subplot(1,3,3)
        pcolor(X, Y, Zeta2)
        axis([0 50 0 50]);
        
        M(istep) = getframe;
    end %if
end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);

