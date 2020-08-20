%== 获得初始时间，设置数据显示格式:
time0 = clock();
format long;

%-- 设置单元尺寸:

Nx = 300;
Ny = 300;
NxNy = Nx * Ny;

dx = 0.03;
dy = 0.03;

%--- 设置积分步:

nstep = 4000;
nprint = 50;
dtime = 1.0e-4;

%--- 材料相关属性:

tau = 0.0003; %时间演化系数
epsilonb = 0.01; %各向异性梯度能量系数平均值0.01
kappa = 1.8; %无量纲潜热，正比于潜热，反比于冷却强度  1.8
delta = 0.02; %各向异性强度0.02
aniso = 6.0; %j：各向异性模数4&6
alpha = 0.9; %过冷系数0.9
gamma = 10.0; %温差放大系数10.0
teq = 1.0; %平衡温度1.0
theta0 = 0.2; %初始偏移角0.2
seed = 10.0; %定义晶核种子大小5.0

%--相/温度初始化:
%[phi, tempr] = nucleus(Nx, Ny, seed);
phi = zeros(Nx, Ny);
tempr = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        if ((i - Nx / 2)^2 + (j - Ny / 2)^2 < seed)
            phi(i, j) = 1.0;
        end
    end
end

%--准备结果存储空间:

X = 0:0.03:9 - 0.03;
Y = 0:0.03:9 - 0.03;
Z = zeros(300, 300);

%- 迭代

for istep = 1:nstep

    %----
    % calculate the laplacians and epsilon:
    %---

    for i = 1:Nx
        for j = 1:Ny
            jp = j + 1;
            jm = j - 1;

            ip = i + 1;
            im = i - 1;

            if (im == 0)
                im = Nx;
            end

            if (ip == (Nx + 1))
                ip = 1;
            end

            if (jm == 0)
                jm = Ny;
            end

            if (jp == (Ny + 1))
                jp = 1;
            end

            hne = phi(ip, j);
            hnw = phi(im, j);
            hns = phi(i, jm);
            hnn = phi(i, jp);
            hnc = phi(i, j);

            lap_phi(i, j) = (hnw + hne + hns + hnn -4.0 * hnc) / (dx * dy);

            hne = tempr(ip, j);
            hnw = tempr(im, j);
            hns = tempr(i, jm);
            hnn = tempr(i, jp);
            hnc = tempr(i, j);

            lap_tempr(i, j) = (hnw + hne + hns + hnn -4.0 * hnc) / (dx * dy);

            %--phidx，phidy:

            phidx(i, j) = (phi(ip, j) - phi(im, j)) / (2.0 * dx);
            phidy(i, j) = (phi(i, jp) - phi(i, jm)) / (2.0 * dy);

            %-- 计算偏移角,求角度用atan2,求值用atan

            theta = atan2(phidy(i, j), phidx(i, j));

            %--- epsl及其导数:

            epsilon(i, j) = epsilonb * (1.0 + delta * cos(aniso * (theta - theta0)));

            epsilon_deriv(i, j) = -epsilonb * aniso * delta * sin(aniso * (theta - theta0));
            %depsl/dtheta

        end %for j

    end %for i

    %----求解相/温度

    for i = 1:Nx

        for j = 1:Ny

            jp = j + 1;
            jm = j - 1;

            ip = i + 1;
            im = i - 1;

            if (im == 0)
                im = Nx;
            end

            if (ip == (Nx + 1))
                ip = 1;
            end

            if (jm == 0)
                jm = Ny;
            end

            if (jp == (Ny + 1))
                jp = 1;
            end

            phiold = phi(i, j);

            %-- first term:

            term1 = (epsilon(i, jp) * epsilon_deriv(i, jp) * phidx(i, jp) - ...
                epsilon(i, jm) * epsilon_deriv(i, jm) * phidx(i, jm)) / (2.0 * dy);

            %-- second term:

            term2 = -(epsilon(ip, j) * epsilon_deriv(ip, j) * phidy(ip, j) - ...
                epsilon(im, j) * epsilon_deriv(im, j) * phidy(im, j)) / (2.0 * dx);

            %-- factor m:

            %m = alpha/pi * atan(gamma*(teq-tempr(i,j)));
            m = alpha / pi * atan(gamma * (teq - tempr(i, j)));
            %-- Time integration:

            phi(i, j) = phi(i, j) +(dtime / tau) * (term1 + term2 + epsilon(i, j)^2 * lap_phi(i, j) + ...
                phiold * (1.0 - phiold) * (phiold -0.5 + m));

            if (phi(i, j) < 10^ - 30)
                phi(i, j) = 10^ - 30;
            end

            if (phi(i, j) > 10^30)
                phi(i, j) = 10^30;
            end

            Z(i, j) = phi(i, j);
            %-- evolve temperature:

            tempr(i, j) = tempr(i, j) +dtime * lap_tempr(i, j) + kappa * (phi(i, j) - phiold);

            if (tempr(i, j) < 10^ - 30)
                tempr(i, j) = 10^ - 30;
            end

            if (tempr(i, j) > 10^30)
                tempr(i, j) = 10^30;
            end

        end

    end

    %---- 输出结果

    if (mod(istep, nprint) == 0)
        subplot(1,2,1)
        pcolor(X, Y, Zcon)
        subplot(1,2,2)
        pcolor(X, Y, Zeta)
        M(istep) = getframe;
    end %if

end %istep

%--- 显示结束时间:

compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);
