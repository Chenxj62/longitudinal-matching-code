function C = peak_slice_sc_coeff_round(infile)

 Nz = 32;
 Nx = 32;
 Ny = 32;

eps0 = 8.854187817e-12;

%% 1) 读取文件
[np] = textread(infile,'%f%*s%*s%*s',1,'headerlines',3);

[x, xp, y, yp, z, p] = textread(infile, ...
    '%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f', ...
    np, 'headerlines', 9); %#ok<ASGLU>

x = x(:); y = y(:); z = z(:);

%% 2) 找峰值电流切片
[counts, edges] = histcounts(z, Nz);
[~, iz] = max(counts);

sel = (z >= edges(iz)) & (z < edges(iz+1));

if nnz(sel) < max(30,20)
    zc = 0.5 * (edges(iz) + edges(iz+1));
    dz = edges(2) - edges(1);
    sel = abs(z - zc) <= 1.5 * dz;
end

xs = x(sel);
ys = y(sel);

if numel(xs) < 30
    error('峰值切片粒子数太少，建议减小 Nz 或适当放宽切片范围。');
end

%% 3) 去质心，并缩放到 sigx = sigy = 1
xs = xs - mean(xs);
ys = ys - mean(ys);

sx = std(xs, 1);
sy = std(ys, 1);

if sx <= 0 || sy <= 0
    error('切片的 sigma_x 或 sigma_y 为 0。');
end

xs = xs / sx;
ys = ys / sy;

%% 4) 建 mesh
xmax = max(6, 1.2 * prctile(abs(xs), 99.5));
ymax = max(6, 1.2 * prctile(abs(ys), 99.5));

xedges = linspace(-xmax, xmax, Nx+1);
yedges = linspace(-ymax, ymax, Ny+1);

dx = xedges(2) - xedges(1);
dy = yedges(2) - yedges(1);

Nxy = histcounts2(xs, ys, xedges, yedges);

% 轻微平滑，减弱 shot noise
kernel = [1 2 1; 2 4 2; 1 2 1] / 16;
Nxy = conv2(Nxy, kernel, 'same');

% 归一化成单位线电荷密度：\int rho dxdy = 1
rho = Nxy / sum(Nxy(:)) / (dx * dy);

xc = 0.5 * (xedges(1:end-1) + xedges(2:end));
yc = 0.5 * (yedges(1:end-1) + yedges(2:end));
[X, Y] = ndgrid(xc, yc); %#ok<NASGU>

%% 5) FFT 解 2D Poisson 求 Ex
Nx2 = 2 * Nx;
Ny2 = 2 * Ny;

rho_pad = zeros(Nx2, Ny2);
ix = floor(Nx/2) + (1:Nx);
iy = floor(Ny/2) + (1:Ny);
rho_pad(ix, iy) = rho;

Lx2 = Nx2 * dx;
Ly2 = Ny2 * dy;

kx = (2*pi/Lx2) * [0:(Nx2/2-1), -Nx2/2:-1];
ky = (2*pi/Ly2) * [0:(Ny2/2-1), -Ny2/2:-1];
[KX, KY] = ndgrid(kx, ky);

rho_hat = fft2(rho_pad);
k2 = KX.^2 + KY.^2;

phi_hat = zeros(size(rho_hat));
mask = (k2 ~= 0);
phi_hat(mask) = rho_hat(mask) ./ (eps0 * k2(mask));

Ex_hat = zeros(size(rho_hat));
Ex_hat(mask) = -1i * KX(mask) .* phi_hat(mask);

Ex_pad = real(ifft2(Ex_hat));
Ex = Ex_pad(ix, iy);

%% 6) 只算系数 C
% 这里 rho 已按单位线电荷归一化，所以 lambda_s = 1
xEx = sum(sum(rho .* X .* Ex)) * dx * dy;

% 圆高斯基准下，C = 8*pi*eps0*<xEx>，且高斯 C = 1
C = 8 * pi * eps0 * xEx;

fprintf('SC shape coefficient (round-scaled, Gaussian = 1): %.10f\n', C);

end