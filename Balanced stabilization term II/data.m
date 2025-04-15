clc,clear
%% 问题描述
% 求解线性抛物方程 u_t = -m_0 * [ - lambda * ( u_xx + u_yy ) + lambda/eta^2 * g(u) ]  + f；周期边界
% 时间上：BDF2-ESAV method
% 空间上：谱配置方法

%% 模型参数
% 求解区域
T = 50; % 时间上
T0 = 0; % 初始时刻
ax = 0;  bx = 1; % 空间上
ay = 0;  by = 1;

epsion = 0.01; % 界面宽度
S = 8.02; % 稳定化参数
C0 = 0; % SAV常数

%% 时间网格剖分
Nt = 5 * T; % 时间步数
tau = (T-T0)/Nt; % 时间步长
t = (T0+tau:tau:T)'; % 时间节点

%% 空间网格剖分
Nx = 2^7;  Ny = Nx; % 谱点数
hx = (bx-ax)/Nx;  hy = (by-ay)/Ny; % 网格长度
x = (ax:hx:bx-hx)';  y = (ay:hy:by-hy)'; % 网格节点

%% 问题初值
[Y,X] = meshgrid(x,x);

% u_2 = -0.8 + ( 0.8 - (-0.8) ) * rand( Nx, Ny );

% save u_intal_06.mat u_2;

uu = load('u_intal_08.mat');
u_2 = uu.u_2;

% uu = load('u_intal_256.mat');
% u_2 = uu.u_2;

u = reshape( u_2, Nx * Ny, 1 );

figure(1);
[X, Y] = meshgrid(y, x);
contourf(X, Y, u_2, 5)
colormap( 'jet' );
colorbar;
axis off

% 记录最大的解u
u_max(1,1) = max( abs(u) );
