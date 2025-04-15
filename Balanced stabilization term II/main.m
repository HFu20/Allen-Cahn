%% 主程序
% 作者：张秉印
% 日期：2023/07/12
clc,clear

%% 导入数据
data;

t1 = clock;
%% 生成刚度矩阵
e = ones( Nx, 1 );
Ax = (epsion^2)/(hx^2) * spdiags([e -2*e e], -1:1, Nx, Ny);
Ax(1,end) = (epsion^2)/(hx^2);
Ax(end,1) = (epsion^2)/(hx^2);
Ix = spdiags(e, 0, Nx, Ny);
A = kron( Ax, Ix ) + kron( Ix, Ax );
I = kron( Ix, Ix );
% full(A)
xi(1,1) = 1;
V(1,1) = 1;

%% 初始能量
[ G, dG] = Fun_Diff( u );
Ef(1,1) = TrapezFun( hx, hy, G ) + C0;
r(1,1) =  Ef(1,1); % 非线性势能
E_inner(1,1) = -1/2 * TrapezFun( hx, hy, (A * u).*u ); % 内能
E(1,1) = E_inner(1,1) + Ef(1,1);

%% 当 n = 1 时，向后欧拉sESAV格式
% 预估步——一阶IMEX格式
% 左侧矩阵
A_left = ( 1 + S * tau ) * I - tau * A;
% 右端项
F = u + tau * S * u + tau * dG;
% 解方程组
u_hat = A_left\F;
% 计算xi
[ G_hat, dG_hat] = Fun_Diff( u_hat );
Ef_hat = TrapezFun( hx, hy, G_hat ) + C0;
xi(2,1) = exp( r(1,1) )/exp( Ef_hat );

% 校正步——一阶sESAV格式
% 左侧矩阵
V(2,1) = FunV( xi(2,1) );
A_left = ( 1 + S * tau ) * I - tau * A;
% 右端项
F = u + tau * S * u_hat + tau * V(2,1) * dG_hat;
% 解u——线性方程组
v = u;
u = A_left\F;
% 解r——线性方程
r(2,1) = r(1,1) - TrapezFun( hx, hy, ( V(2,1) * dG_hat - S * ( u - u_hat ) ).*( u - v ) ); 

% 记录最大的解u
u_max(2,1) = max( abs(u) );

%% n >= 2时，BDF2-sESAV格式
for n = 2:Nt
    % 能量
    [ G, dG] = Fun_Diff( u );

    % 离散的能量
    Ef(n,1) = TrapezFun( hx, hy, G ) + C0; % 非线性势能
    E_inner(n,1) = -1/2 * TrapezFun( hx, hy, (A * u).*u ); % 内能
    E(n,1) = E_inner(n,1) + Ef(n,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 预估步——一阶IMEX格式
    % 左侧矩阵
    A_left = ( 1 + S * tau ) * I - tau * A;

    % 右端项
    F = u + tau * S * u + tau * dG;
    % 解方程组
    u_hat = A_left\F;
%     u = u_hat;
    % 计算xi
    [ G_hat, dG_hat] = Fun_Diff( u_hat );
    Ef_hat = TrapezFun( hx, hy, G_hat ) + C0;
    xi(n+1,1) = exp( r(n,1) )/exp( Ef_hat );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 校正步——BDF2-sESAV格式
    % 左侧矩阵
    V(n+1,1) = FunV( xi(n+1,1) );
    A_left = ( 3 + S * 2 * tau ) * I - 2 * tau * A;
    % 右端项
    F = 4 * u - v + S * 2 * tau * u_hat + 2 * tau * V(n+1,1) * dG_hat;
    % 解方程组
    v = u;
    u = A_left\F;
    % 解r——线性方程
    r(n+1,1) = r(n,1) - TrapezFun( hx, hy, ( V(n+1,1) * dG_hat - S * ( u - u_hat ) ).*( u - v ) );

    % 记录最大的解u
    u_max(n+1,1) = max( abs(u) );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% 离散的能量
Ef(n+1,1) = TrapezFun( hx, hy, G ) + C0; % 非线性势能
E_inner(n+1,1) = -1/2 * TrapezFun( hx, hy, (A * u).*u ); % 内能
E(n+1,1) = E_inner(n+1,1) + Ef(n+1,1);

%% 终止时间
t2 = clock;
UPC_times = etime(t2, t1);

% save( 'u_N8.mat', 'u' );

%% 数值解
figure(2);
u_2 = reshape(u, Nx, Ny);
% u_2(end+1,:) = u_2(1,:);
% u_2(:,end+1) = u_2(:,1);
contourf(X, Y, u_2, 5)
colormap( 'jet' );
colorbar;
%axis off

%% 能量耗散律
figure(3);
tt(1,1) = 0;
tt(2:Nt+1,1) = t;
plot( tt, E, 'k-.', 'LineWidth', 1 )

%% 解的最大值
figure(4);
plot( [0,T]', [0.9575,0.9575]', 'r--', 'LineWidth', 1 )
hold on
plot( tt, u_max, 'k-.', 'LineWidth', 1 )

% plot( tt, u_max, 'k-.', 'LineWidth', 1 )
% axis( [0 50 1 1+2*10^(-3)] )

%% xi
figure(5);
plot( [0,T]', [1,1]', 'r--', 'LineWidth', 1 )
hold on
plot( tt, xi, 'k-.', 'LineWidth', 1 )

figure(6);
plot( [0,T]', [1,1]', 'r--', 'LineWidth', 1 )
hold on
plot( tt, V, 'k-.', 'LineWidth', 1 )

%% 数据存储
save E_11_tau02.mat E;
save uT50_11_tau02.mat u_2;
save u_max_11_tau02.mat u_max;

% save xi_11_01.mat xi;
% save V_11_01.mat V;


% save( 'u_T1.mat', 'u' );
