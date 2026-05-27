% RUN_ACC_ARC_SSDR_MP  Example: ss + dr CSR through an arc
%
%   Same pipeline as run_acc_arc_mp.m, but uses the combined
%   steady-state + drift-CSR feature of csr_calc_mp:
%     - one call with 'drift_csr' = true   -> ss + dr (total)
%     - one call with 'drift_csr' = false  -> ss only (for comparison)
%   Benchmarked against BMAD particle dumps in paper4/.

clear; clc; close all;

addpath('..\..');       % readBmad, getElementMatrix, createElement

% =================== USER INPUTS ===================
bunchfile   = '..\DOG3_ENT_200K_V13a.ast';            % particle dump at injection
Q_total     = 70e-12;                    % [C]

% --- Cavity (one-pass acceleration) -----------------------------------
E_ref0   = 1000;                          % [MeV] in
E_target = 1000;                          % [MeV] out
phi_deg  = 13.55;                         % [deg] off-crest
f        = 1.3e9;                         % [Hz]
V_cav    = (E_target - E_ref0) / cosd(phi_deg);

% --- Arc ---
latticefile = '..\arc3wc.bmad';
T566        = [];                       % [m]

% --- CSR ---
n_grid         = 200;
smooth_win     = 40;
n_substeps     = 100;
sigzf_override = [];                      % [] = auto

% --- BMAD benchmark ---
benchfile_nc1 = '..\..\paper4\pap4outnc.dat';
benchfile_nc = '..\..\paper4\pap4outwc.dat';
benchfile_wc = '..\..\paper4\pap4outwc2.dat';
n_bins_out   = 100;

font_name = 'Cambria';
% ===================================================

% --- 1) load particles -----------------------------------------------
[z0, delta0] = read_dist(bunchfile);
N_part = numel(z0);
fprintf('Loaded %d particles\n', N_part);

% --- 2) cavity -------------------------------------------------------
[z1, delta1, E_ref1] = cavity_long(z0, delta0, E_ref0, V_cav, phi_deg, f);


% --- 3) lattice ------------------------------------------------------
beamline  = readBmad(latticefile);
gamma_rel = E_ref1 / 0.510998950;
r=nonlinear_calculator(beamline,5,6,6,0.0001)
r=calculateTotalMatrix(beamline);
r(5,6)

% --- 4) ss + dr CSR (combined, the new default) ---------------------
out_total = csr_calc_mp(z1, delta1, Q_total, beamline, gamma_rel, T566, ...
        'n_grid',     n_grid, ...
        'smooth_win', smooth_win, ...
        'n_substeps', n_substeps, ...
        'sigzf',      sigzf_override, ...
        'drift_csr',  true, ...
        'plot_flag',  1);

% --- 5) ss only (for comparison against the combined result) --------
out_ss = csr_calc_mp(z1, delta1, Q_total, beamline, gamma_rel, T566, ...
        'n_grid',     n_grid, ...
        'smooth_win', smooth_win, ...
        'n_substeps', n_substeps, ...
        'sigzf',      sigzf_override, ...
        'drift_csr',  false, ...
        'plot_flag',  0);

% --- 6) BMAD benchmark beams ----------------------------------------
c           = 2.99792458e8;
q_per_macro = Q_total / N_part;



[zc_tot, I_tot] = current_from_particles(out_total.z_out_wCSR,  n_bins_out, q_per_macro, c);

% --- 7) plots (subplot, no titles, Cambria) -------------------------
figure('Color','w','Position',[100 100 700 300]);

% phase space
subplot(1,2,1);
% plot(zb_nc1*1e6, db_nc1*1e3, 'k.', 'MarkerSize',2);hold on;
% plot(zb_nc*1e6, db_nc*1e3, '+', 'Color',[.55 .75 1.0], 'MarkerSize',2); 
% plot(zb_wc*1e6, db_wc*1e3, '+', 'Color',[1.0 .55 .55], 'MarkerSize',2);
% plot(out_ss.z_out_wCSR   *1e6, out_ss.delta_wCSR   *1e3, 'b.', 'MarkerSize',2);
plot(out_total.z_out_wCSR*1e6, out_total.delta_wCSR*1e3,'.','Color',[.8 .0 .0], 'MarkerSize',2);
grid on;
xlabel('z  [\mum]'); ylabel('\delta  [10^{-3}]');
legend({'BMAD no CSR','BMAD ss-CSR','BMAD all CSR','analysis ss','analysis ss+dr'}, ...
       'Location','best');
set(gca,'FontName',font_name); ax = gca; ax.LineWidth = 1.5;

% current
subplot(1,2,2);
% plot(zc_no *1e6, I_no,  'k--',  'LineWidth',1.4);hold on;
% plot(zc_bnc*1e6, I_bnc, 'b--', 'LineWidth',1.2); 
% plot(zc_bwc*1e6, I_bwc, 'Color',[.8 .0 .0],'LineStyle','--', 'LineWidth',1.2);
% plot(zc_ss *1e6, I_ss,  'b-',  'LineWidth',1.4);
plot(zc_tot*1e6, I_tot, 'Color',[.8 .0 .0],  'LineWidth',1.4);
grid off;
xlabel('z  [\mum]'); ylabel('I  [A]');
legend({'BMAD no CSR','BMAD ss-CSR','BMAD all CSR','analysis ss','analysis ss+dr'}, ...
       'Location','best');
set(gca,'FontName',font_name); ax = gca; ax.LineWidth = 1.5;




% =========================================================================
% Local helpers
% =========================================================================
function [z, delta] = read_dist(infile)
    data = dlmread(infile);
    if size(data, 2) >= 10
        data = data(data(:, 10) > 0, :);
    end
    if size(data, 1) < 2
        error('read_dist: %s has < 2 valid rows.', infile);
    end
    z_ref = data(1, 3);
    p_ref = data(1, 6);
    z_abs = data(2:end, 3) + z_ref;
    p_abs = data(2:end, 6) + p_ref;
    keep  = isfinite(z_abs) & isfinite(p_abs);
    z     = z_abs(keep);
    p     = p_abs(keep);
    z     = z - mean(z);
    pbar  = mean(p);
    delta = (p - pbar) / pbar;
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

function [z, delta] = read_bmad_zdelta(infile)
    raw = readmatrix(infile, 'FileType','text', ...
                     'NumHeaderLines', 9, 'ConsecutiveDelimitersRule','join');
    if size(raw,2) < 6
        error('read_bmad_zdelta: %s has only %d columns.', infile, size(raw,2));
    end
    keep  = isfinite(raw(:,5)) & isfinite(raw(:,6));
    z     = raw(keep, 5);
    delta = raw(keep, 6);
end
