clear; clc; close all;

%% =========================================================
% User input
%% =========================================================
infile = 'fvv12.beam';

%% =========================================================
% Options
%% =========================================================
Nz_select   = 20;      % find peak-current slice
Ngrid       = 16;      % 2D mesh (IGF is accurate even at low Ngrid)
eta_min     = 0.4;
eta_max     = 6.0;
Neta        = 120;
minSliceNum = 80;

margin_uv   = 1.02;    % mesh box from FULL bunch edge in normalized coords

%% =========================================================
% 1) Read file
%% =========================================================
fprintf('Reading file: %s\n', infile);

[np] = textread(infile,'%f%*s%*s%*s',1,'headerlines',3);
[nc] = textread(infile,'%f%*s%*s%*s',1,'headerlines',6); %#ok<NASGU>
[x, xp, y, yp, z, p] = textread(infile, ...
    '%f%f%f%f%f%f%*f%*f%*f%*f%*f%*f%*f', ...
    np, 'headerlines', 9); %#ok<ASGLU>

x = x(:);
y = y(:);
z = z(:);

fprintf('Total particles read = %d\n', numel(x));

%% =========================================================
% 2) Find peak-current slice
%% =========================================================
[counts, edges] = histcounts(z, Nz_select);
[~, iz] = max(counts);

sel = (z >= edges(iz)) & (z < edges(iz+1));

if nnz(sel) < max(120, round(0.002*np))
    zc = 0.5 * (edges(iz) + edges(iz+1));
    dzs = edges(2) - edges(1);
    sel = abs(z - zc) <= 1.5 * dzs;
end

if nnz(sel) < minSliceNum
    error('Too few particles in selected peak slice.');
end

xs_raw = x(sel);
ys_raw = y(sel);

fprintf('Particles in peak slice = %d\n', numel(xs_raw));

%% =========================================================
% 3) Use peak-slice centroid as origin
%% =========================================================
x0 = mean(xs_raw);
y0 = mean(ys_raw);

x_all = x - x0;
y_all = y - y0;

xs = xs_raw - x0;
ys = ys_raw - y0;

sigx0 = std(xs, 1);
sigy0 = std(ys, 1);
eta0  = sigy0 / sigx0;

if sigx0 <= 0 || sigy0 <= 0
    error('Invalid sigma_x or sigma_y.');
end

fprintf('Peak-slice sigma_x0 = %.10e\n', sigx0);
fprintf('Peak-slice sigma_y0 = %.10e\n', sigy0);
fprintf('Physical eta0       = %.10f\n', eta0);

%% =========================================================
% 4) Normalize using PEAK-SLICE sigmas
%    u,v for selected slice
%    u_all,v_all for FULL bunch edge definition
%% =========================================================
u = xs / sigx0;
v = ys / sigy0;

u_all = x_all / sigx0;
v_all = y_all / sigy0;

fprintf('rms(u_slice)=%.8f, rms(v_slice)=%.8f\n', sqrt(mean(u.^2)), sqrt(mean(v.^2)));

%% =========================================================
% 5) Mesh box from FULL bunch edge in normalized coords
%% =========================================================
umax = margin_uv * max(abs(u_all));
vmax = margin_uv * max(abs(v_all));

u_edges = linspace(-umax, umax, Ngrid+1);
v_edges = linspace(-vmax, vmax, Ngrid+1);

du = u_edges(2) - u_edges(1);
dv = v_edges(2) - v_edges(1);

fprintf('\nMesh defined by FULL bunch edge (normalized):\n');
fprintf('  Ngrid = %d\n', Ngrid);
fprintf('  umax  = %.10e\n', umax);
fprintf('  vmax  = %.10e\n', vmax);
fprintf('  du    = %.10e\n', du);
fprintf('  dv    = %.10e\n', dv);

%% =========================================================
% 6) Eta scan — IGF (Integrated Green Function) solver
%
% Uses the OpenSC method:
%   - CIC deposit
%   - 2D Integrated Green Functions for direct field calc
%   - FFT convolution (zero-padded for open boundary)
%   - CIC gather
%
% Self-normalization: C = X_pic + Y_pic removes all constant factors
%% =========================================================
eta_list = linspace(eta_min, eta_max, Neta);

X_pic  = zeros(size(eta_list));
Y_pic  = zeros(size(eta_list));
Mx_pic = zeros(size(eta_list));
My_pic = zeros(size(eta_list));
X_th   = zeros(size(eta_list));
Y_th   = zeros(size(eta_list));
Rx     = zeros(size(eta_list));
Ry     = zeros(size(eta_list));

for m = 1:Neta
    eta = eta_list(m);

    out = solve_igf_2d(u, v, eta, u_edges, v_edges);

    % Field moments
    X_pic(m) = mean(u  .* out.Eu_p);
    Y_pic(m) = mean(-v .* out.dphidv_p);

    % Gaussian theory (fraction)
    X_th(m) = 1 / (1 + eta);
    Y_th(m) = eta / (1 + eta);

    % Self-normalize
    C = X_pic(m) + Y_pic(m);
    Mx = X_pic(m) / C;
    My = Y_pic(m) / C;

    Mx_pic(m) = Mx;
    My_pic(m) = My;

    Rx(m) = Mx * (1 + eta);
    Ry(m) = My * (1 + eta) / eta;

    fprintf('eta = %.4f | X_pic = %.8e | Y_pic = %.8e | Rx = %.8f | Ry = %.8f\n', ...
        eta, X_pic(m), Y_pic(m), Rx(m), Ry(m));
end

%% =========================================================
% 7) Interpolate at actual physical eta0
%% =========================================================
X_pic_0 = interp1(eta_list, X_pic, eta0, 'pchip', 'extrap');
Y_pic_0 = interp1(eta_list, Y_pic, eta0, 'pchip', 'extrap');

X_th_0  = 1 / (1 + eta0);
Y_th_0  = eta0 / (1 + eta0);

C0  = X_pic_0 + Y_pic_0;
Rx0 = X_pic_0 / C0 * (1 + eta0);
Ry0 = Y_pic_0 / C0 * (1 + eta0) / eta0;

fprintf('\n====================================================\n');
fprintf('Eta-scan result at physical eta0\n');
fprintf('eta0   = %.10f\n', eta0);
fprintf('X_pic  = %.12e\n', X_pic_0);
fprintf('Y_pic  = %.12e\n', Y_pic_0);
fprintf('X_th   = %.12e\n', X_th_0);
fprintf('Y_th   = %.12e\n', Y_th_0);
fprintf('Rx0    = %.12f\n', Rx0);
fprintf('Ry0    = %.12f\n', Ry0);
fprintf('====================================================\n');

%% =========================================================
% 8) Diagnostics at eta0
%% =========================================================
diag_out = solve_igf_2d(u, v, eta0, u_edges, v_edges);

%% =========================================================
% 9) Plot
%% =========================================================
figure('Color','w','Position',[70 40 1600 320]);

subplot(2,3,1);
histogram(u, 120, 'Normalization', 'pdf');
hold on;
uu = linspace(-8,8,600);
plot(uu, exp(-0.5*uu.^2)/sqrt(2*pi), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('u = x/\sigma_{x0}');
ylabel('PDF');
title('Peak slice u projection');

subplot(2,3,2);
histogram(v, 120, 'Normalization', 'pdf');
hold on;
plot(uu, exp(-0.5*uu.^2)/sqrt(2*pi), 'r--', 'LineWidth', 1.5);
grid on;
xlabel('v = y/\sigma_{y0}');
ylabel('PDF');
title('Peak slice v projection');

subplot(2,3,3);
imagesc(diag_out.u_c, diag_out.v_c, diag_out.rho.');
axis xy image;
xlabel('u');
ylabel('v');
title('Peak-slice density on FULL-bunch-defined mesh');
colorbar;

subplot(1,3,1);
plot(eta_list, Mx_pic, 'b-', 'LineWidth', 1.8); hold on;
plot(eta_list, X_th,   'k--', 'LineWidth', 1.5);
plot(eta0, Rx0/(1+eta0), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(eta0, X_th_0,       'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('\eta');
ylabel('M_x(\eta)');
%title('M_x = X_{pic}/C  vs  1/(1+\eta)');
legend('IGF ','Gaussian theory','Location','best');

subplot(1,3,2);
plot(eta_list, My_pic, 'm-', 'LineWidth', 1.8); hold on;
plot(eta_list, Y_th,   'k--', 'LineWidth', 1.5);
plot(eta0, Ry0*eta0/(1+eta0), 'mo', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(eta0, Y_th_0,            'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('\eta');
ylabel('M_y(\eta)');
%title('M_y = Y_{pic}/C  vs  \eta/(1+\eta)');
legend('IGF ','Gaussian theory','Location','best');

subplot(1,3,3);
plot(eta_list, Rx, 'b-', 'LineWidth', 1.8); hold on;
plot(eta_list, Ry, 'm--', 'LineWidth', 1.8);
yline(1, 'k:', 'LineWidth', 1.2);
plot(eta0, Rx0, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(eta0, Ry0, 'mo', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('\eta');
ylabel('R(\eta) = IGF / Gaussian');
%title('R_x, R_y  (= 1 for Gaussian)');
legend('R_x(\eta)','R_y(\eta)','Location','best');

sgtitle(sprintf('IGF solver | file: %s | Ngrid=%d', infile, Ngrid), ...
    'Interpreter','none');

%% =========================================================
% 10) Export
%% =========================================================
result = struct();
result.infile = infile;
result.n_total = numel(x);
result.n_slice = numel(u);

result.sigx0 = sigx0;
result.sigy0 = sigy0;
result.eta0  = eta0;

result.Ngrid = Ngrid;
result.umax  = umax;
result.vmax  = vmax;
result.du    = du;
result.dv    = dv;

result.eta_list = eta_list;
result.X_pic  = X_pic;
result.Y_pic  = Y_pic;
result.Mx_pic = Mx_pic;
result.My_pic = My_pic;
result.X_th   = X_th;
result.Y_th   = Y_th;
result.Rx     = Rx;
result.Ry     = Ry;

result.X_pic_0 = X_pic_0;
result.Y_pic_0 = Y_pic_0;
result.X_th_0  = X_th_0;
result.Y_th_0  = Y_th_0;
result.Rx0 = Rx0;
result.Ry0 = Ry0;

fprintf('\nDone.\n');

%% =========================================================
% Local functions
%% =========================================================

function out = solve_igf_2d(u, v, eta, u_edges, v_edges)
% 2D Integrated Green Function solver (OpenSC method)
%
% Direct field calculation: E_u, E_y computed from IGF field kernels
% via FFT convolution.  No finite-difference of phi needed.
%
% IGF advantage over point-source Green's function:
%   - Cell-integrated kernel is exact (no self-term hack)
%   - Much more accurate at low Ngrid
%   - Self-term automatically zero for field kernels (by symmetry)

    u = u(:);
    v = v(:);

    Nu = numel(u_edges) - 1;
    Nv = numel(v_edges) - 1;

    du = u_edges(2) - u_edges(1);
    dv = v_edges(2) - v_edges(1);

    u_c = 0.5 * (u_edges(1:end-1) + u_edges(2:end));
    v_c = 0.5 * (v_edges(1:end-1) + v_edges(2:end));

    % ---------------------------
    % CIC deposit
    % ---------------------------
    rho = cic_deposit_2d(u, v, u_edges, v_edges);
    rho = rho / (sum(rho(:)) * du * dv);

    % ---------------------------
    % Zero-padded arrays for open-boundary convolution
    % ---------------------------
    Npu = 2 * Nu;
    Npv = 2 * Nv;

    rho_pad = zeros(Npu, Npv);
    rho_pad(1:Nu, 1:Nv) = rho;

    % Wrapped grid shifts (cell-center separations)
    iu = 0:Npu-1;
    iv = 0:Npv-1;
    diu = iu;
    div = iv;
    diu(iu >= Nu) = iu(iu >= Nu) - Npu;
    div(iv >= Nv) = iv(iv >= Nv) - Npv;

    [II, JJ] = ndgrid(diu, div);
    DU = II * du;       % u-separation between cell centers
    DV = JJ * dv;       % v-separation between cell centers

    % ---------------------------
    % Cell boundaries in physical coordinates
    %   physical x = u,  physical y = eta * v
    % ---------------------------
    x1 = DU - du/2;          x2 = DU + du/2;
    y1 = eta*(DV - dv/2);    y2 = eta*(DV + dv/2);

    cell_area = du * eta * dv;   % physical cell area

    % ---------------------------
    % IGF for E_u = -dphi/du  (2D analogue of OpenSC xlafun)
    %
    % Kernel: x / (2*pi*(x^2 + y^2))
    % Antiderivative: xlafun2d(x,y) = y*ln(r) - y + x*atan(y/x)
    %
    % Self-term (DU=0,DV=0) is automatically 0 (antisymmetric kernel)
    % ---------------------------
    Gex = (xlafun2d(x2,y2) - xlafun2d(x1,y2) ...
         - xlafun2d(x2,y1) + xlafun2d(x1,y1)) / (2*pi * cell_area);

    % ---------------------------
    % IGF for E_y = -dphi/dy  (2D analogue of OpenSC ylafun)
    %
    % Kernel: y / (2*pi*(x^2 + y^2))
    % Antiderivative: ylafun2d(x,y) = x*ln(r) - x + y*atan(x/y)
    % ---------------------------
    Gey = (ylafun2d(x2,y2) - ylafun2d(x1,y2) ...
         - ylafun2d(x2,y1) + ylafun2d(x1,y1)) / (2*pi * cell_area);

    % ---------------------------
    % FFT convolution  (field = kernel * rho * du * dv)
    % ---------------------------
    rho_fft = fft2(rho_pad);
    Eu_pad = real(ifft2(rho_fft .* fft2(Gex))) * du * dv;
    Ey_pad = real(ifft2(rho_fft .* fft2(Gey))) * du * dv;

    Eu = Eu_pad(1:Nu, 1:Nv);          % E_u = -dphi/du
    Ey = Ey_pad(1:Nu, 1:Nv);          % E_y = -dphi/dy (physical y)

    % dphi/dv = dphi/dy * dy/dv = (-Ey) * eta
    dphidv = -Ey * eta;

    % ---------------------------
    % CIC gather back to particles
    % ---------------------------
    Eu_p     = cic_gather_2d(Eu,     u, v, u_edges, v_edges);
    dphidv_p = cic_gather_2d(dphidv, u, v, u_edges, v_edges);

    out = struct();
    out.rho      = rho;
    out.Eu       = Eu;
    out.dphidv   = dphidv;
    out.Eu_p     = Eu_p;
    out.dphidv_p = dphidv_p;
    out.u_c      = u_c;
    out.v_c      = v_c;
end

% ---------------------------------------------------------------
% 2D antiderivative of x/(x^2+y^2) — analogous to OpenSC xlafun
%
%   xlafun2d(x,y) = y * ln(sqrt(x^2+y^2)) - y + x * atan(y/x)
%
% Satisfies: d^2/dxdy [xlafun2d] = x/(x^2+y^2)
% ---------------------------------------------------------------
function val = xlafun2d(x, y)
    r2 = x.^2 + y.^2;

    % ln(r): safe at r=0
    lnr = zeros(size(r2));
    nz = r2 > 0;
    lnr(nz) = 0.5 * log(r2(nz));

    % x*atan(y/x): limit is 0 when x->0
    xatanyx = zeros(size(x));
    mx = (x ~= 0);
    xatanyx(mx) = x(mx) .* atan(y(mx) ./ x(mx));

    val = y .* lnr - y + xatanyx;
    val(~nz) = 0;   % (0,0) -> 0
end

% ---------------------------------------------------------------
% 2D antiderivative of y/(x^2+y^2) — analogous to OpenSC ylafun
%
%   ylafun2d(x,y) = x * ln(sqrt(x^2+y^2)) - x + y * atan(x/y)
%
% Satisfies: d^2/dxdy [ylafun2d] = y/(x^2+y^2)
% ---------------------------------------------------------------
function val = ylafun2d(x, y)
    r2 = x.^2 + y.^2;

    lnr = zeros(size(r2));
    nz = r2 > 0;
    lnr(nz) = 0.5 * log(r2(nz));

    yatanxy = zeros(size(y));
    my = (y ~= 0);
    yatanxy(my) = y(my) .* atan(x(my) ./ y(my));

    val = x .* lnr - x + yatanxy;
    val(~nz) = 0;
end

% ---------------------------------------------------------------
% CIC (Cloud-In-Cell) deposit: bilinear weight to 4 neighboring cells
% ---------------------------------------------------------------
function rho = cic_deposit_2d(u, v, u_edges, v_edges)

    u = u(:);
    v = v(:);
    Np = numel(u);

    Nu = numel(u_edges) - 1;
    Nv = numel(v_edges) - 1;

    du = u_edges(2) - u_edges(1);
    dv = v_edges(2) - v_edges(1);

    u_c = 0.5 * (u_edges(1:end-1) + u_edges(2:end));
    v_c = 0.5 * (v_edges(1:end-1) + v_edges(2:end));

    fu = (u - u_c(1)) / du + 1;
    fv = (v - v_c(1)) / dv + 1;

    iu0 = floor(fu);
    iv0 = floor(fv);

    wu1 = fu - iu0;
    wv1 = fv - iv0;
    wu0 = 1 - wu1;
    wv0 = 1 - wv1;

    iu0 = max(1, min(Nu, iu0));
    iu1 = max(1, min(Nu, iu0 + 1));
    iv0 = max(1, min(Nv, iv0));
    iv1 = max(1, min(Nv, iv0 + 1));

    rho = zeros(Nu, Nv);
    for p = 1:Np
        rho(iu0(p), iv0(p)) = rho(iu0(p), iv0(p)) + wu0(p) * wv0(p);
        rho(iu1(p), iv0(p)) = rho(iu1(p), iv0(p)) + wu1(p) * wv0(p);
        rho(iu0(p), iv1(p)) = rho(iu0(p), iv1(p)) + wu0(p) * wv1(p);
        rho(iu1(p), iv1(p)) = rho(iu1(p), iv1(p)) + wu1(p) * wv1(p);
    end
end

% ---------------------------------------------------------------
% CIC (Cloud-In-Cell) gather: bilinear interpolation from 4 cells
% ---------------------------------------------------------------
function fp = cic_gather_2d(F, u, v, u_edges, v_edges)

    u = u(:);
    v = v(:);

    Nu = numel(u_edges) - 1;
    Nv = numel(v_edges) - 1;

    du = u_edges(2) - u_edges(1);
    dv = v_edges(2) - v_edges(1);

    u_c = 0.5 * (u_edges(1:end-1) + u_edges(2:end));
    v_c = 0.5 * (v_edges(1:end-1) + v_edges(2:end));

    fu = (u - u_c(1)) / du + 1;
    fv = (v - v_c(1)) / dv + 1;

    iu0 = floor(fu);
    iv0 = floor(fv);

    wu1 = fu - iu0;
    wv1 = fv - iv0;
    wu0 = 1 - wu1;
    wv0 = 1 - wv1;

    iu0 = max(1, min(Nu, iu0));
    iu1 = max(1, min(Nu, iu0 + 1));
    iv0 = max(1, min(Nv, iv0));
    iv1 = max(1, min(Nv, iv0 + 1));

    ind00 = sub2ind([Nu, Nv], iu0, iv0);
    ind10 = sub2ind([Nu, Nv], iu1, iv0);
    ind01 = sub2ind([Nu, Nv], iu0, iv1);
    ind11 = sub2ind([Nu, Nv], iu1, iv1);

    fp = wu0.*wv0.*F(ind00) + wu1.*wv0.*F(ind10) + ...
         wu0.*wv1.*F(ind01) + wu1.*wv1.*F(ind11);
end
