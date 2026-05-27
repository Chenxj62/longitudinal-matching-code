% RUN_ACC_ARC  Example: one-pass cavity 10 MeV -> 1000 MeV, then arc CSR
%
%   Pipeline:
%     1. bmad2slice(bunchfile)         -> (z, delta, I)  at 10 MeV
%     2. cavity_long(... 1.3 GHz ...)  -> (z, delta)     at 1000 MeV
%     3. readBmad('../../pap4/latfv2.arc') -> beamline
%     4. csr_calc_slice(...)           -> end-of-arc phase space + I_out

clear; clc; close all;

addpath('..');           % bmad2slice.m   (in CSRwake/)
addpath('..\..');        % readBmad.m, getElementMatrix.m, createElement.m

% =================== USER INPUTS ===================
%bunchfile   = '..\..\paper4\pap4.beam';%'..\beam_merger.ast';      % particle dump at 10 MeV stage
bunchfile   = '..\DOG3_ENT_200K_V13a.ast';
Q_total     = 70e-12;                    % [C]

% --- Cavity (one-pass acceleration) -----------------------------------
E_ref0   = 1000;                           % [MeV]   in
E_target = 1000;                         % [MeV]   out
phi_deg  = 13.55;                            % [deg]   off-crest (0 = on-crest)
f        = 1.3e9;                        % [Hz]    RF frequency

% required voltage to reach E_target on this phase
V_cav    = (E_target - E_ref0) / cosd(phi_deg);

% --- Arc (lattice file) -----------------------------------------------
latticefile = '..\..\paper4\latfv2.arc';
latticefile = '..\fvv13.arc';
T566        =-1.37;                         % [m]  2nd-order longitudinal dispersion

% --- Slice extraction (bmad2slice) ------------------------------------
nslices       = 300;
fit_order     = 9;
min_particles = 10;

% --- CSR settings -----------------------------------------------------
smooth_win     = 30;
n_substeps     = 30;
n_bends        = 0;                      % 0 = use all detected
dR56           = -0.0008;
sigzf_override = [];                     % [] = auto-derive from chirp

% --- BMAD benchmark ---------------------------------------------------
benchfile_nc  = '..\..\paper4\pap4outnc.dat';   % BMAD: no CSR
benchfile_ss  = '..\..\paper4\pap4outwc.dat';   % BMAD: ss CSR only
benchfile_all = '..\..\paper4\pap4outwc2.dat';  % BMAD: ss + dr (all) CSR
n_bins_out    = 80;                            % bins for benchmark current

font_name = 'Cambria';
% ===================================================

% --- 1) slice arrays at 10 MeV ----------------------------------------
[z0, delta0, I0] = bmad2slice(bunchfile, nslices, ...
                              fit_order, min_particles, Q_total);



% --- 2) cavity: 10 MeV -> 1000 MeV ------------------------------------
[z1, delta1, E_ref1] = cavity_long(z0, delta0, E_ref0, V_cav, phi_deg, f);



% --- 3) arc lattice ---------------------------------------------------
beamline  = readBmad(latticefile);
r=calculateTotalMatrix(beamline);
r=r(5,6)+dR56;
fprintf('r56=%f\n',r);
gamma_rel = E_ref1 / 0.510998950;

% --- 4a) ss + dr CSR (combined, the default) --------------------------
out_total = csr_calc_slice(z1, delta1, I0, beamline, gamma_rel, T566, ...
        'smooth_win', smooth_win, ...
        'n_substeps', n_substeps, ...
        'n_bends',    n_bends, ...
        'dR56',       dR56, ...
        'sigzf',      sigzf_override, ...
        'drift_csr',  true, ...
        'plot_flag',  1);

% --- 4b) ss only (for comparison) -------------------------------------
out_ss = csr_calc_slice(z1, delta1, I0, beamline, gamma_rel, T566, ...
        'smooth_win', smooth_win, ...
        'n_substeps', n_substeps, ...
        'n_bends',    n_bends, ...
        'dR56',       dR56, ...
        'sigzf',      sigzf_override, ...
        'drift_csr',  false, ...
        'plot_flag',  0);

% --- 5) BMAD benchmark beams and currents -----------------------------
c = 2.99792458e8;
% [zb_nc,  db_nc ] = read_bmad_zdelta(benchfile_nc);
% [zb_ss,  db_ss ] = read_bmad_zdelta(benchfile_ss);
% [zb_all, db_all] = read_bmad_zdelta(benchfile_all);
% qpm_b_nc  = 1e-10 / numel(zb_nc);
% qpm_b_ss  = 1e-10 / numel(zb_ss);
% qpm_b_all = 1e-10 / numel(zb_all);
% 
% [zc_bnc,  I_bnc ] = current_from_particles(zb_nc,  n_bins_out, qpm_b_nc,  c);
% [zc_bss,  I_bss ] = current_from_particles(zb_ss,  n_bins_out, qpm_b_ss,  c);
% [zc_ball, I_ball] = current_from_particles(zb_all, n_bins_out, qpm_b_all, c);

% --- 5b) slice output currents via Jacobian I = I0 / |dz_out/dz_0| ----
zp_ss   = gradient(out_ss.z_out_wCSR,     out_ss.z0);
zp_tot  = gradient(out_total.z_out_wCSR,  out_total.z0);
I_ss    = I0 ./ max(abs(zp_ss),  eps);
I_tot   = I0 ./ max(abs(zp_tot), eps);

% --- 6) plots (subplot, no titles, Cambria) ---------------------------
figure('Color','w','Position',[100 100 700 300]);

% phase space
subplot(1,2,1);
%plot(zb_nc *1e6, db_nc *1e3, '+', 'Color','k', 'MarkerSize',2); hold on;
%plot(zb_ss *1e6, db_ss *1e3, '+', 'Color',[.55 .75 1.0], 'MarkerSize',2);
%plot(zb_all*1e6, db_all*1e3, '+', 'Color',[1.0 .55 .55], 'MarkerSize',2);
%plot(out_ss.z_out_wCSR    *1e6, out_ss.delta_wCSR    *1e3, 'b.-', ...
%     'LineWidth',1.0, 'MarkerSize',6);
plot(out_total.z_out_wCSR *1e6, out_total.delta_wCSR *1e3, ...
     'Color',[.8 .0 .0], 'LineStyle','-', 'Marker','.', ...
     'LineWidth',1.0, 'MarkerSize',6);
grid off;
xlabel('z  [\mum]'); ylabel('\delta  [10^{-3}]');
% legend({'Bmad no CSR','Bmad ss CSR','Bmad all CSR', ...
%         'analysis ss','analysis ss+dr'}, ...
%        'Location','best');
set(gca,'FontName',font_name); ax = gca; ax.LineWidth = 1.5;

% current
subplot(1,2,2);
%plot(zc_bnc *1e6, I_bnc,  'k+', 'LineWidth',1.2); hold on;
%plot(zc_bss *1e6, I_bss,  'b+', 'LineWidth',1.2);
%plot(zc_ball*1e6, I_ball, 'r+', 'LineWidth',1.2);
%plot(out_ss.z_out_wCSR   *1e6, I_ss,  'b-',  'LineWidth',1.4);
plot(out_ss.z_out_wCSR*1e6, I_tot, 'Color',[.8 .0 .0], ...
     'LineStyle','-', 'LineWidth',1.4);
grid off;
xlabel('z  [\mum]'); ylabel('I  [A]');
% legend({'Bmad no CSR','Bmad ss CSR','Bmad all CSR', ...
%         'analysis ss','analysis ss+dr'}, ...
%        'Location','best');
set(gca,'FontName',font_name); ax = gca; ax.LineWidth = 1.5;




% =========================================================================
% Local helpers
% =========================================================================
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
% Read z (col 5) and delta (col 6) from a BMAD ASCII beam file.
% Skip 9 header lines, drop any non-numeric rows.
    raw = readmatrix(infile, 'FileType','text', ...
                     'NumHeaderLines', 9, 'ConsecutiveDelimitersRule','join');
    if size(raw,2) < 6
        error('read_bmad_zdelta: %s has only %d columns.', infile, size(raw,2));
    end
    keep  = isfinite(raw(:,5)) & isfinite(raw(:,6));
    z     = raw(keep, 5);
    delta = raw(keep, 6);
end
