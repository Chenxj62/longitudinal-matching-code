% RUN_CSR_MP  Example driver: particle beam + R56 line -> CSR transport.
%
%   Multi-particle counterpart of run_csr.m. The output current is built
%   directly by histogramming the transported particle distribution
%   z_out_wCSR, NOT from the slice Jacobian I0 / |dz_out/dz_0|. The slice
%   approach breaks down whenever dz_out/dz_0 changes sign (folding) and
%   is noisy near caustics; the particle histogram does not care.

 clear; clc; close all;

% =================== USER INPUTS ===================
bunchfile = 'Bunch.ast';          % ASTRA-style particle file
Q_total   = 70e-12;               % [C] total bunch charge
r56file   = 'r56v13.line_dat';    % R56(s) line file (cols: idx, s, R56)

gamma_rel = 1000 / 0.511;         % 1 GeV electrons -> gamma ~ 1957

% bend radii in s-order (10 bends, same as run_csr.m)
rho_list = [ 4.3*ones(1,4), 2.79*ones(1,4), -3.8, 3.8 ];

% sigma_zf override (auto-derived if empty). sige is auto = std(delta).
sigzf = 12e-6;

% second-order longitudinal dispersion of the line
T566 = -1.7;

% lambda(z) grid + Gaussian smoothing inside csr_calc_mp
n_grid     = 300;
smooth_win = 40;

% number of bends to integrate over (0 = all detected)
n_bends = 0;

% dR56 added ONLY to the transport step (no effect on the line integrals)
dR56 = 0.002;

% Optional symmetric sigma-cut to drop tail particles before transport
% (set to 0 to keep all particles).
sigma_cut = 0;

% Histogram resolution for the output-current diagnostic
n_bins_out = 400;
% ===================================================

% --- 1) load particles ----------------------------------------------------
[z0, delta0] = load_particles_ast(bunchfile);

if sigma_cut > 0
    sz   = std(z0);
    keep = abs(z0 - mean(z0)) < sigma_cut * sz;
    z0     = z0(keep);
    delta0 = delta0(keep);
end

N_part = numel(z0);
fprintf('Loaded %d particles, Q_total = %.3g pC, sigma_z = %.3g um.\n', ...
        N_part, Q_total*1e12, std(z0)*1e6);

% --- 2) end-to-end CSR + transport (per particle) ------------------------
out = csr_calc_mp(z0, delta0, Q_total, r56file, gamma_rel, rho_list, ...
        T566, ...
        'n_grid',     n_grid, ...
        'smooth_win', smooth_win, ...
        'n_bends',    n_bends, ...
        'dR56',       dR56, ...
        's_col',      2, ...
        'R56_col',    3, ...
        'plot_flag',  1, ...
        'sigzf',      sigzf);

% --- 3) output current by histogramming transported particles -----------
% q_per_macro * N_bin * c / dz   (same convention as bmad2slice.m)
c           = 2.99792458e8;
q_per_macro = Q_total / N_part;

[zc_in,  I_in]  = current_from_particles(z0,              n_bins_out, q_per_macro, c);
[zc_no,  I_no]  = current_from_particles(out.z_out_noCSR, n_bins_out, q_per_macro, c);
[zc_w,   I_w ]  = current_from_particles(out.z_out_wCSR,  n_bins_out, q_per_macro, c);

figure('Name','run\_csr\_mp current','Color','w','Position',[100 100 1200 750]);

subplot(2,2,1);
plot(zc_in*1e6, I_in, 'k-', 'LineWidth',1.2); grid on;
xlabel('z  [\mum]'); ylabel('I(z)  [A]');
title(sprintf('Input current  (N=%d particles, %d bins)', N_part, n_bins_out));

subplot(2,2,2);
plot(zc_no*1e6, I_no, 'b-', 'LineWidth',1.4); hold on;
plot(zc_w *1e6, I_w,  'r-', 'LineWidth',1.4); grid on;
xlabel('z_{out}  [\mum]'); ylabel('I_{out}  [A]');
legend('no CSR','with CSR','Location','best');
title('Output current (histogram of transported particles)');

subplot(2,2,3);
plot(out.z_out_noCSR*1e6, out.delta_noCSR*1e3, '.', ...
     'Color',[.55 .75 1.0], 'MarkerSize',2); hold on;
plot(out.z_out_wCSR *1e6, out.delta_wCSR *1e3, '.', ...
     'Color',[1.0 .45 .45], 'MarkerSize',2);
grid on;
xlabel('z_{out}  [\mum]'); ylabel('\delta_{out}  [10^{-3}]');
legend('no CSR','with CSR','Location','best');
title('End-of-line phase space');

% CSR kicks are not stored in out; compute on the fly.
dz_csr     = out.z_out_wCSR  - out.z_out_noCSR;
ddelta_csr = out.delta_wCSR  - out.delta_noCSR;

subplot(2,2,4);
plot(z0*1e6, dz_csr*1e6,     'k.', 'MarkerSize',2); hold on;
plot(z0*1e6, ddelta_csr*1e3, 'r.', 'MarkerSize',2);
grid on;
xlabel('z_0  [\mum]');
legend('dz_{csr}  [\mum]','d\delta_{csr}  [10^{-3}]','Location','best');
title('Per-particle CSR kicks');

ud = out.Properties.UserData;
fprintf('\nResolved scalars:\n');
fprintf('  sigma_z0 = %.3g um   sigma_zf = %.3g um\n', ud.sigma_z0*1e6, ud.sigma_zf*1e6);
fprintf('  sige_used = %.3g     r56_tot = %.3g m\n',   ud.sige_used, ud.r56_tot);
fprintf('  S_R56 = %.3g  S_one = %.3g\n',              ud.S_R56, ud.S_one);
fprintf('\nPeaks  I_in = %.2f A   I_out(no CSR) = %.2f A   I_out(w/ CSR) = %.2f A\n', ...
        max(I_in), max(I_no), max(I_w));


% =========================================================================
% Local helpers
% =========================================================================
function [z, delta] = load_particles_ast(infile)
% LOAD_PARTICLES_AST  Read an ASTRA-style particle dump.
%   Mirrors the ASTRA branch of bmad2slice.m:
%     - dlmread the file
%     - filter rows where column 10 > 0 (status flag)
%     - row 1 is the reference particle (absolute values)
%     - rows 2..end store z and momentum as deltas relative to row 1
%   Returns z [m] (centered) and delta = (p - <p>) / <p>.

    data = dlmread(infile);
    data = data(data(:, 10) > 0, :);
    if size(data, 1) < 2
        error('load_particles_ast: %s has < 2 valid rows.', infile);
    end

    z_ref = data(1, 3);
    p_ref = data(1, 6);

    z_abs = data(2:end, 3) + z_ref;     % rows 2..end are relative to row 1
    p_abs = data(2:end, 6) + p_ref;

    z    = z_abs - mean(z_abs);
    pbar = mean(p_abs);
    delta = (p_abs - pbar) / pbar;
end

function [zc, I] = current_from_particles(z, nb, q_per_macro, c)
    zmin  = min(z);
    zmax  = max(z);
    edges = linspace(zmin, zmax, nb + 1);
    dz    = edges(2) - edges(1);
    zc    = 0.5 * (edges(1:end-1) + edges(2:end));
    cnt   = histcounts(z, edges);
    I     = q_per_macro * cnt * c / dz;
    zc    = zc(:);
    I     = I(:);
end
