function result = MBI_calc(beamline, beam, opts)
% MBI_calc  Microbunching instability gain — CSR + optics (SI units)
%
% Self-contained Volterra integral solver.
% Based on C.-Y. Tsai, PRAB 19, 114401 (2016).
%
% ALL UNITS ARE SI (meters, radians, seconds, amperes).
%
% Inputs:
%   beamline  - struct array from createElement
%   beam      - struct:
%       .energy_MeV                          [MeV]
%       .I_b                                 [A]
%       .emit_norm_x, .emit_norm_y           [m-rad]
%       .betax0, .betay0                     [m]
%       .alphax0, .alphay0                   [-]
%       .sigma_delta                         [-]
%       .chirp (optional)                    [m^-1]
%       .sigma_z                             [m]
%   opts      - struct (all optional):
%       .lambda_range     [m_start, m_end]   (default: auto)
%       .scan_num         (default: 50)
%       .mesh_num         (default: 400)
%       .iCSR_ss          steady-state CSR   (default: 1)
%       .iCSR_tr          transient CSR      (default: 0)
%       .iCSR_drift       CSR drift effect   (default: 0)
%       .iLSC             LSC impedance      (default: 0)
%       .LSC_model        1=on-axis, 2=avg, 3=Gaussian (default: 3)
%       .issCSRpp         CSR shielding      (default: 0)
%       .full_pipe_height pipe height [m]    (default: 1e50)
%       .itransLD         transverse LD      (default: 1)
%       .calc_energy_gain                    (default: 1)
%       .plot                                (default: 0)

%% Defaults
if nargin < 3, opts = struct(); end
def = struct('lambda_range',[], 'scan_num',50, 'mesh_num',400, ...
    'iCSR_ss',1, 'iCSR_tr',0, 'iCSR_drift',0, ...
    'iLSC',0, 'LSC_model',3, ...
    'issCSRpp',0, 'full_pipe_height',1e50, ...
    'itransLD',1, 'calc_energy_gain',1, 'plot',0);
for f = fieldnames(def)'
    if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
end
if ~isfield(beam,'chirp'), beam.chirp = 0; end

%% Constants (SI)
me_MeV = 0.511;
re     = 2.81794e-15;       % classical electron radius [m]
I_A    = 17045;              % Alfven current [A]
A_csr  = 3^(-1/3) * gamma(2/3) * (1i*sqrt(3) - 1);

gamma0 = beam.energy_MeV / me_MeV;
ex     = beam.emit_norm_x / gamma0;   % geometric emittance [m-rad]
ey     = beam.emit_norm_y / gamma0;
bx0    = beam.betax0;                  % [m]
by0    = beam.betay0;
ax0    = beam.alphax0;
ay0    = beam.alphay0;
sd     = beam.sigma_delta;
chirp  = beam.chirp;                   % [m^-1]
sig_z  = beam.sigma_z;                 % [m]
nb     = beam.I_b / (re * I_A);

% ======================================================================
%  Step 1: Walk beamline — cumulative R(0->s), gamma(s), dipole table
% ======================================================================
nelem = length(beamline);
s_arr = zeros(nelem+1, 1);            % position [m]
R_all = cell(nelem+1, 1);
g_arr = zeros(nelem+1, 1);

R_all{1} = eye(6);  g_arr(1) = gamma0;
R_cum = eye(6);  cur_g = gamma0;
dip_ent = []; dip_ext = []; dip_rho = [];

for ie = 1:nelem
    [Re, go] = getElementMatrix(beamline(ie), cur_g);
    R_cum = Re * R_cum;  cur_g = go;
    s_arr(ie+1) = s_arr(ie) + beamline(ie).length;
    R_all{ie+1} = R_cum;
    g_arr(ie+1) = cur_g;
    if strcmpi(beamline(ie).type, 'dipole') && abs(beamline(ie).h) > 1e-10
        dip_ent(end+1) = ie;                   %#ok
        dip_ext(end+1) = ie + 1;               %#ok
        dip_rho(end+1) = 1 / beamline(ie).h;   %#ok  [m]
    end
end
Npts = nelem + 1;

% ======================================================================
%  Step 2: Extract R5x(0->s) — all in SI [m]
% ======================================================================
R51 = zeros(Npts,1); R52 = R51; R53 = R51;
R54 = R51; R55 = R51; R56 = R51;
for ip = 1:Npts
    R = R_all{ip};
    R51(ip) = R(5,1); R52(ip) = R(5,2);
    R53(ip) = R(5,3); R54(ip) = R(5,4);
    R55(ip) = R(5,5); R56(ip) = R(5,6);
end

% HK sign convention
R51 = -R51;  R52 = -R52;  R53 = -R53;  R54 = -R54;  R56 = -R56;
% Velocity bunching R56 (disabled — not needed for magnetic chicane)
% R56 = R56 + (s_arr - s_arr(1)) ./ g_arr.^2;
% Compression factor (beam-based, as in Tsai's volterra_mat code)
% C = sigma_z(0) / sqrt(R56^2*sd^2 + (R55-h*R56)^2*sigma_z(0)^2)
% This prevents C from diverging near full compression;
% energy spread limits the achievable compression.
C = sig_z ./ sqrt(R56.^2 * sd^2 + (R55 - chirp * R56).^2 * sig_z^2);

% Beam sizes along beamline (from Twiss propagation)
% sigma11 = beta*eps, sigma12 = -alpha*eps, sigma22 = gamma_tw*eps
% where gamma_tw = (1+alpha^2)/beta
gx0 = (1 + ax0^2) / bx0;
gy0 = (1 + ay0^2) / by0;
sig_x_arr = zeros(Npts, 1);
sig_y_arr = zeros(Npts, 1);
for ip = 1:Npts
    Rx = R_all{ip}(1:2, 1:2);
    Ry = R_all{ip}(3:4, 3:4);
    % beta_x(s) from Twiss transport: sigma_11 = R11^2*beta0 - 2*R11*R12*alpha0 + R12^2*gamma0
    bx_s = Rx(1,1)^2*bx0 - 2*Rx(1,1)*Rx(1,2)*ax0 + Rx(1,2)^2*gx0;
    by_s = Ry(1,1)^2*by0 - 2*Ry(1,1)*Ry(1,2)*ay0 + Ry(1,2)^2*gy0;
    sig_x_arr(ip) = sqrt(ex * bx_s);
    sig_y_arr(ip) = sqrt(ey * by_s);
end

% ======================================================================
%  Step 3: Dipole table [m]
% ======================================================================
N_D = length(dip_ent);
dip = struct('N', N_D, ...
    's1', zeros(N_D,1), 's2', zeros(N_D,1), 'rho', zeros(N_D,1));
for id = 1:N_D
    dip.s1(id)  = s_arr(dip_ent(id));
    dip.s2(id)  = s_arr(dip_ext(id));
    dip.rho(id) = dip_rho(id);          % [m]
end

% ======================================================================
%  Step 4: Mesh and precompute
% ======================================================================
N  = opts.mesh_num;
s0 = s_arr(1);  s1 = s_arr(end) - 1e-3;   % small offset [m]
sm = linspace(s0, s1, N)';
hs = (s1 - s0) / N;                        % mesh step [m]

% Remove duplicate s points (from zero-length elements like markers/edges)
[s_u, iu] = unique(s_arr, 'stable');
C_u = C(iu); R51u = R51(iu); R52u = R52(iu);
R53u = R53(iu); R54u = R54(iu); R55u = R55(iu); R56u = R56(iu);
sx_u = sig_x_arr(iu); sy_u = sig_y_arr(iu); g_u = g_arr(iu);

Cm   = interp1(s_u, C_u,   sm, 'linear', 'extrap');
R51m = interp1(s_u, R51u, sm, 'linear', 'extrap');
R52m = interp1(s_u, R52u, sm, 'linear', 'extrap');
R53m = interp1(s_u, R53u, sm, 'linear', 'extrap');
R54m = interp1(s_u, R54u, sm, 'linear', 'extrap');
R55m = interp1(s_u, R55u, sm, 'linear', 'extrap');
R56m = interp1(s_u, R56u, sm, 'linear', 'extrap');
sx_m = interp1(s_u, sx_u, sm, 'linear', 'extrap');
sy_m = interp1(s_u, sy_u, sm, 'linear', 'extrap');
gm   = interp1(s_u, g_u,  sm, 'linear', 'extrap');

% --- Debug plot: effective current along beamline ---
figure; plot(sm, Cm * beam.I_b, 'r-', 'LineWidth', 2);
xlabel('s (m)'); ylabel('I_{eff} = I_b \cdot C(s)  [A]'); grid on;
title('Effective bunch current along beamline');
fprintf('DEBUG: I_b = %.2f A, C at exit = %.4f, I_eff at exit = %.2f A\n', ...
    beam.I_b, Cm(end), Cm(end)*beam.I_b);

rho_m = ones(N,1) * 1e50;
for id = 1:N_D
    mask = (sm >= dip.s1(id)) & (sm <= dip.s2(id));
    rho_m(mask) = dip.rho(id);
end

% ======================================================================
%  Step 5: Auto wavelength range [m]
% ======================================================================
if isempty(opts.lambda_range)
    R56end = abs(R56m(end));
    if R56end < 1e-6
        lopt = 3 * mean(abs(R56m)) * sd;
    else
        lopt = 3 * R56end * sd;
    end
    lam_start = max(0.2*lopt, 1e-7);
    lam_end   = max(50*lopt, 1e-2);
else
    lam_start = opts.lambda_range(1);
    lam_end   = opts.lambda_range(2);
end
Nscan  = opts.scan_num;
iEcalc = opts.calc_energy_gain;

% ======================================================================
%  Step 6: Volterra solver — wavelength scan
% ======================================================================
lam_arr = linspace(lam_start, lam_end, Nscan);
Gf_dd  = zeros(1, Nscan);
G_vs_s = zeros(N, Nscan);
if iEcalc
    Gf_de = zeros(1,Nscan); Gf_ed = zeros(1,Nscan);
    Gf_ee = zeros(1,Nscan);
end

coeff0 = 1i * re * nb / gamma0;

fprintf('MBI_calc: %d wavelengths [%.4e, %.4e] m ...\n', ...
    Nscan, lam_start, lam_end);

for iscan = 1:Nscan
    lam = lam_arr(iscan);              % [m]
    k = 2*pi / lam;                    % [m^-1]

    % --- Landau damping at mesh ---
    LDx = 0.5*k^2*ex/bx0 * Cm.^2 .* ...
          (bx0^2*(R51m - R52m*ax0/bx0).^2 + R52m.^2);
    LDy = 0.5*k^2*ey/by0 * Cm.^2 .* ...
          (by0^2*(R53m - R54m*ay0/by0).^2 + R54m.^2);
    LDd = 0.5*k^2*sd^2 * Cm.^2 .* R56m.^2;
    if opts.itransLD, LD = LDx+LDy+LDd; else, LD = LDd; end

    % --- Initial gain vectors ---
    n_1k0 = 1.0;  e_1k0 = 0.0;
    G0_k  = n_1k0 * exp(-LD);
    G0_kp = -1i*k * Cm .* R56m * e_1k0 .* exp(-LD);
    E0_k  = -1i*k * Cm .* R56m * sd^2 * n_1k0 .* exp(-LD);
    E0_kp = (1 - k^2*Cm.^2.*R56m.^2*sd^2) * e_1k0 .* exp(-LD);

    % --- Build kernel matrix ---
    K_mat  = zeros(N);
    ML_mat = zeros(N);

    for j = 1:(N-1)
        Cj = Cm(j);

        % Total impedance Z at s'
        Z = compute_Z_csr(k, Cj, sm(j), rho_m(j), dip, gamma0, opts, A_csr);
        if opts.iLSC
            Z = Z + lsc_impedance(k*Cj, sx_m(j), sy_m(j), gm(j), opts.LSC_model);
        end
        if Z == 0, continue; end

        ii  = (j+1):N;
        fac = 1;  if j == 1, fac = 0.5; end

        % Modified transport dR5x(s,s') = C(s)*R5x(s) - C(s')*R5x(s')
        dR51 = Cm(ii).*R51m(ii) - Cj*R51m(j);
        dR52 = Cm(ii).*R52m(ii) - Cj*R52m(j);
        dR53 = Cm(ii).*R53m(ii) - Cj*R53m(j);
        dR54 = Cm(ii).*R54m(ii) - Cj*R54m(j);
        dR56 = Cm(ii).*R56m(ii) - Cj*R56m(j);

        % R56(s'->s) symplectic
        R56t = R56m(ii) - R56m(j)*R55m(ii)/R55m(j) ...
             + R51m(j)*R52m(ii) - R51m(ii)*R52m(j) ...
             + R53m(j)*R54m(ii) - R53m(ii)*R54m(j);

        % Landau damping in kernel
        LDx_k = 0.5*k^2*ex/bx0 * (bx0^2*(dR51-dR52*ax0/bx0).^2 + dR52.^2);
        LDy_k = 0.5*k^2*ey/by0 * (by0^2*(dR53-dR54*ay0/by0).^2 + dR54.^2);
        LDd_k = 0.5*k^2*sd^2 * dR56.^2;
        if opts.itransLD, LD_k = LDx_k+LDy_k+LDd_k; else, LD_k = LDd_k; end

        % Density-density kernel
        K_mat(ii, j) = fac * coeff0*k .* Cm(ii) .* Cj .* Z ...
                      .* R56t .* exp(-LD_k);

        if iEcalc
            K_M = fac * re*nb/gamma0 * k^2 .* Cm(ii).^2 .* Cj ...
                * sd^2 .* Z .* dR56 .* R56t .* exp(-LD_k);
            K_L = fac * re*nb/gamma0 .* Cj .* Z .* exp(-LD_k);
            ML_mat(ii, j) = K_M - K_L;
        end
    end

    % --- Solve Volterra: (I - h*K) \ G0 ---
    I_K = eye(N) - hs * K_mat;
    gkd = I_K \ G0_k;
    G_c_dd = gkd / max(abs(G0_k(1)), 1e-30);
    Gf_dd(iscan)     = abs(G_c_dd(end));
    G_vs_s(:, iscan) = abs(G_c_dd);

    if iEcalc
        gkp = I_K \ G0_kp;
        ekd = hs * ML_mat * (I_K \ G0_k)  + E0_k;
        ekp = hs * ML_mat * (I_K \ G0_kp) + E0_kp;
        Gf_ed(iscan) = abs(gkp(end) / max(abs(e_1k0), 1e-30));
        Gf_de(iscan) = abs(ekd(end) / max(abs(G0_k(1)), 1e-30));
        Gf_ee(iscan) = abs(ekp(end) / max(abs(e_1k0), 1e-30));
    end

    if mod(iscan, max(floor(Nscan/10),1)) == 0 || iscan == Nscan
        [pk_s, pk_idx] = max(abs(G_c_dd));
        fprintf('  [%d/%d] lam=%.2f um, Gf_end=%.4e, Gf_max=%.4e at s=%.2f m\n', ...
            iscan, Nscan, lam_arr(iscan)*1e6, Gf_dd(iscan), pk_s, sm(pk_idx));
    end
end

% ======================================================================
%  Output
% ======================================================================
result.lambda  = lam_arr;          % [m]
result.s       = sm;               % [m]
result.Gf_dd   = Gf_dd;
result.G_vs_s  = G_vs_s;
% Peak gain at any position along beamline (not just exit)
[result.Gf_dd_max, max_s_idx] = max(G_vs_s, [], 1);
result.s_max = sm(max_s_idx);     % position of peak gain for each lambda
if iEcalc
    result.Gf_de    = Gf_de;
    result.Gf_ed    = Gf_ed;
    result.Gf_ee    = Gf_ee;
    result.Gf_total = abs(Gf_dd + Gf_ed);
    result.Ef_total = abs(Gf_de + Gf_ee);
end

% ======================================================================
%  Plot (display lambda axis in um for readability)
% ======================================================================
if opts.plot && Nscan >= 2
    lam_um = lam_arr * 1e6;
    figure; plot(lam_um, Gf_dd, 'b-', 'LineWidth', 2);
    xlabel('\lambda (\mum)'); ylabel('Density gain'); grid on;
    title('MBI: density-to-density');
    if iEcalc
        figure; plot(lam_um, result.Gf_total, 'k-', 'LineWidth', 2);
        xlabel('\lambda (\mum)'); ylabel('Total density gain'); grid on;
        title('MBI: total density gain');
    end
    figure; surf(lam_um, result.s, G_vs_s);
    xlabel('\lambda (\mum)'); ylabel('s (m)'); zlabel('|G|');
    title('MBI gain map'); shading interp; colorbar;
end

[pk, pidx] = max(Gf_dd);
fprintf('MBI_calc done.\n');
fprintf('  Peak gain at EXIT:     %.4e at lambda = %.2f um\n', pk, lam_arr(pidx)*1e6);
[pk2, pidx2] = max(result.Gf_dd_max);
fprintf('  Peak gain ANYWHERE:    %.4e at lambda = %.2f um, s = %.2f m\n', ...
    pk2, lam_arr(pidx2)*1e6, result.s_max(pidx2));
end
