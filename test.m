infile='fvv12.beam';
Nz = 16; 
Nx = 132; 
Ny = 132; 
eps0 = 8.854187817e-12;

%% ---------------------------
% 1) 读取粒子文件
%% ---------------------------
[np] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 3);
[nc] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 6); %#ok<NASGU>
[x1, xp, y1, yp, z, p] = textread(infile, ...
    '%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f', ...
    np, 'headerlines', 9); %#ok<ASGLU>
x = x1(:)/3; y = y1(:); z = z(:);

%% ---------------------------
% 2) 找峰值电流切片（z 上粒子数最多的 bin）
%% ---------------------------
[counts, edges] = histcounts(z, Nz);
[~, iz] = max(counts);
sel = (z >= edges(iz)) & (z < edges(iz+1));

% 如果这个 bin 太窄，自动稍微放宽一点，降低统计噪声

xs = x(sel);
ys = y(sel);
if numel(xs) < 30
    error('峰值切片粒子数太少，建议增大 Nz 或手动加宽 z-window。');
end

%% ---------------------------
% 3) 去掉切片质心
%% ---------------------------
xs = xs - mean(xs);
ys = ys - mean(ys);

%% ---------------------------
% 4) 缩放 xs 使其 RMS 等于 ys
%% ---------------------------
sx_orig = std(xs, 1);
sy_orig = std(ys, 1);
scale_factor = sy_orig / sx_orig;
xs = xs * 1;   % 缩放 x 坐标

%% ---------------------------
% 5) 建二维 mesh
%% ---------------------------
sx0 = std(xs, 1)
sy0 = std(ys, 1)
xmax = max(6*sx0, 1.2*prctile(abs(xs), 99.5));
ymax = max(6*sy0, 1.2*prctile(abs(ys), 99.5));
xedges = linspace(-xmax, xmax, Nx+1);
yedges = linspace(-ymax, ymax, Ny+1);
dx = xedges(2) - xedges(1);
dy = yedges(2) - yedges(1);

% 计数沉积到网格
Nxy = histcounts2(xs, ys, xedges, yedges);
% 做一个轻微平滑，降低 shot noise
kernel = [1 2 1; 2 4 2; 1 2 1] / 16;
Nxy = conv2(Nxy, kernel, 'same');
% 归一化成“单位线电荷”密度：\int rho dxdy = 1
rho = Nxy / sum(Nxy(:)) / (dx*dy);

xc = 0.5 * (xedges(1:end-1) + xedges(2:end));
yc = 0.5 * (yedges(1:end-1) + yedges(2:end));
[X, Y] = ndgrid(xc, yc);

%% ---------------------------
% 6) Hockney method: free-space 2D Green's function convolution求 Ex
%% ---------------------------
Nx2 = 2 * Nx;
Ny2 = 2 * Ny;

% 零填充 rho（放在角落，标准 Hockney 布局）
rho_pad = zeros(Nx2, Ny2);
rho_pad(1:Nx, 1:Ny) = rho;

% 构建 Ex 自由空间 Green 函数: Gx(dx,dy) = dx / (2*pi*eps0*(dx^2+dy^2))
dxi_vec = [(0:Nx-1), (-Nx:-1)] * dx;
dyj_vec = [(0:Ny-1), (-Ny:-1)] * dy;
[DXI, DYJ] = ndgrid(dxi_vec, dyj_vec);
R2 = DXI.^2 + DYJ.^2;
Kx_green = zeros(Nx2, Ny2);
valid = R2 > 0;
Kx_green(valid) = DXI(valid) ./ (2*pi*eps0 * R2(valid));

% 构建 Ey 自由空间 Green 函数: Gy(dx,dy) = dy / (2*pi*eps0*(dx^2+dy^2))
Ky_green = zeros(Nx2, Ny2);
Ky_green(valid) = DYJ(valid) ./ (2*pi*eps0 * R2(valid));

% FFT 卷积 + 离散积分因子 dx*dy
rho_hat = fft2(rho_pad);
Ex_pad = real(ifft2(rho_hat .* fft2(Kx_green))) * dx * dy;
Ey_pad = real(ifft2(rho_hat .* fft2(Ky_green))) * dx * dy;
Ex = Ex_pad(1:Nx, 1:Ny);
Ey = Ey_pad(1:Nx, 1:Ny);

%% ---------------------------
% 7) mesh 上算 <x Ex> 及 rms 宽度
%% ---------------------------
lambda_s = sum(rho(:)) * dx * dy;  % 理论上 = 1
xEx_num = sum(sum(rho .* X .* Ex)) * dx * dy;
yEy_num = sum(sum(rho .* Y .* Ey)) * dx * dy;

sigx = sqrt(sum(sum(rho .* X.^2)) * dx * dy / lambda_s);
sigy = sqrt(sum(sum(rho .* Y.^2)) * dx * dy / lambda_s);

%% ---------------------------
% 8) 高斯参考值及输出结果
%% ---------------------------
xEx_gauss = lambda_s / (4*pi*eps0) * sigx / (sigx + sigy);
yEy_gauss = lambda_s / (4*pi*eps0) * sigy / (sigx + sigy);

ratio_x = xEx_num / xEx_gauss;
ratio_y = yEy_num / yEy_gauss;

fprintf('ratio_x = <xEx>/<xEx>_gauss = %.6f\n', ratio_x);
fprintf('ratio_y = <yEy>/<yEy>_gauss = %.6f\n', ratio_y);