function out = csr_drift_mp(z, delta, Q_total, beamline, gamma_rel, T566, varargin)
% CSR_DRIFT_MP  Drift / transient (exit) CSR kicks for a lattice.
%
%   out = csr_drift_mp(z, delta, Q_total, beamline, gamma_rel, T566, ...
%                      NAME, VALUE, ...)
%
%   Companion to csr_calc_mp.m.  Implements the drift-CSR formulas:
%
%     z_csr_D(u)      = k_csrD * lambda(u) * sum_i R56_i * C_i * log(2 L_Di / (rho phi_i) + 1)
%     delta_csr_D(u)  = k_csrD * lambda(u) * sum_i        C_i * log(2 L_Di / (rho phi_i) + 1)
%
%   with
%     k_csrD  = 2 * r_e * N / (gamma * sigma_zf)
%     phi_0   = (6 * sigma_z_in / |rho|)^(1/3)
%     phi_i   = phi_0 / C_i^(1/3)
%     C_i     = 1 / (1 - (sige/sigzf) * R56_i)     (compression at exit of bend i)
%     R56_i   = downstream-side R56 of the drift segment after bend i
%
%   The downstream-side R56 is obtained the same way as in
%   ../../CalcCSRkick.m: the beamline is reversed and walked forward;
%   R_cum(5,6) at the entry of each drift segment is the residual R56
%   from that point to the end of the (original) line.
%
%   The output has the same shape as csr_calc_mp.m so the two can be
%   combined transparently in a driver.
%
%   Required (same as csr_calc_mp)
%     z, delta, Q_total, beamline, gamma_rel, T566
%
%   Optional name-value pairs (subset of csr_calc_mp)
%     'sige'        ([])     override sige; default = std(delta)
%     'sigzf'       ([])     override sigzf; default = std(z + R56*delta + T566*delta^2)
%     'smooth_win'  (5)      Gaussian smoothing window for binned lambda
%     'n_grid'      (200)    grid points for lambda(z) build
%     'interp'      ('linear') how lambda is interpolated onto particles
%     'dR56'        (0)      added to R56 in the transport step only
%     'plot_flag'   (0)      diagnostic plot
%
%   Output: table with one row per particle, columns
%     z_out_noCSR, z_out_wCSR, delta_noCSR, delta_wCSR
%   where z_out_wCSR / delta_wCSR include ONLY the drift-CSR contribution.
%
%   out.Properties.UserData carries
%     Q_total, N_part, sige_used, R56_file, r56_tot,
%     sigma_z0, sigma_zf, k_csrD, S_drift_one, S_drift_R56,
%     n_segments_used, dz_csr, ddelta_csr (per-particle drift CSR kicks),
%     z_grid, lambda_grid

    % ---------- options ----------
    p = inputParser;
    p.addParameter('sige',       []);
    p.addParameter('sigzf',      []);
    p.addParameter('smooth_win', 5);
    p.addParameter('n_grid',     200);
    p.addParameter('interp',     'linear');
    p.addParameter('dR56',       0);
    p.addParameter('plot_flag',  0);
    p.parse(varargin{:});
    opt = p.Results;

    % ---------- sanitize ----------
    z     = z(:);
    delta = delta(:);
    if numel(z) ~= numel(delta)
        error('csr_drift_mp: z and delta length mismatch.');
    end
    if any(~isfinite(z)) || any(~isfinite(delta))
        error('csr_drift_mp: non-finite z or delta.');
    end
    if ~isstruct(beamline) || isempty(beamline)
        error('csr_drift_mp: beamline must be a non-empty struct array.');
    end
    N_part = numel(z);

    % ---------- moments ----------
    sigz0 = std(z);
    if isempty(opt.sige) || (isscalar(opt.sige) && isnan(opt.sige))
        sige = std(delta);
    else
        sige = opt.sige;
    end

    % ---------- total R56 from the lattice ----------
    R56_file = local_walk_total_R56(beamline);
    r56_tot  = R56_file + opt.dR56;
    if R56_file ~= 0 && opt.dR56 ~= 0
        sige = sige * (R56_file + opt.dR56) / R56_file;
    end

    % ---------- sigzf (no-CSR transport std, with R56*delta + T566*delta^2) ----------
    z_transported = z + r56_tot * delta + T566 * delta.^2;
    sigzf_auto    = std(z_transported);

    if isempty(opt.sigzf) || (isscalar(opt.sigzf) && isnan(opt.sigzf))
        sigzf = sigzf_auto;
        if sigzf <= 0
            error('csr_drift_mp: sigma_zf <= 0.  Past full compression.');
        end
    else
        sigzf = max(opt.sigzf, sigzf_auto);
    end

    % ---------- drift CSR line sums (reversed-walk like CalcCSRkick) ----------
    [S_drift_R56, S_drift_one, n_seg] = local_drift_sums( ...
        beamline, sige, sigzf, sigz0);

    if n_seg == 0
        warning('csr_drift_mp: no drift segments between bends; drift CSR = 0.');
    end

    % ---------- k_csrD prefactor ----------
    r_e    = 2.8179403262e-15;
    qe     = 1.602176634e-19;
    N_elec = Q_total / qe;
    k_csrD = 2 * r_e * N_elec / (gamma_rel * sigzf);

    % ---------- lambda(z) on uniform grid ----------
    [z_grid, lambda_grid] = local_lambda(z, opt.n_grid, opt.smooth_win, sigz0);

    % ---------- interpolate lambda to each particle ----------
    lambda_part = interp1(z_grid, lambda_grid, z, opt.interp, 0);

    % ---------- per-particle drift CSR kicks ----------
    ddelta_csr = k_csrD * lambda_part * S_drift_one;
    dz_csr     = k_csrD * lambda_part * S_drift_R56;

    % ---------- transport (drift CSR only) ----------
    z_out_noCSR = z + r56_tot * delta + T566 * delta.^2;
    z_out_wCSR  = z_out_noCSR + dz_csr;
    delta_noCSR = delta;
    delta_wCSR  = delta + ddelta_csr;

    % ---------- output ----------
    out = table(z_out_noCSR, z_out_wCSR, delta_noCSR, delta_wCSR, ...
                'VariableNames', {'z_out_noCSR','z_out_wCSR', ...
                                  'delta_noCSR','delta_wCSR'});
    out.Properties.UserData = struct( ...
        'Q_total',         Q_total, ...
        'N_part',          N_part, ...
        'sige_used',       sige, ...
        'R56_file',        R56_file, ...
        'r56_tot',         r56_tot, ...
        'sigma_z0',        sigz0, ...
        'sigma_zf',        sigzf, ...
        'k_csrD',          k_csrD, ...
        'S_drift_one',     S_drift_one, ...
        'S_drift_R56',     S_drift_R56, ...
        'n_segments_used', n_seg, ...
        'dz_csr',          dz_csr, ...
        'ddelta_csr',      ddelta_csr, ...
        'z_grid',          z_grid, ...
        'lambda_grid',     lambda_grid );

    % ---------- optional plot ----------
    if opt.plot_flag
        figure('Color','w','Name','csr\_drift\_mp', ...
               'Position',[100 100 1100 450]);

        subplot(1,2,1);
        plot(z_grid*1e6, lambda_grid, 'k-', 'LineWidth',1.4); grid on;
        xlabel('z  [\mum]'); ylabel('\lambda(z)');

        subplot(1,2,2);
        plot(z*1e6, dz_csr*1e6,     'k.', 'MarkerSize',2); hold on;
        plot(z*1e6, ddelta_csr*1e3, 'r.', 'MarkerSize',2); grid on;
        xlabel('z_0  [\mum]');
        legend('dz_{csr,drift}  [\mum]','d\delta_{csr,drift}  [10^{-3}]', ...
               'Location','best');
    end
end


%% =====================================================================
function R56_total = local_walk_total_R56(beamline)
    R_cum = eye(6);
    for k = 1:numel(beamline)
        R_cum = getElementMatrix(beamline(k)) * R_cum;
    end
    R56_total = R_cum(5,6);
end


%% =====================================================================
function [is_bend, rho] = local_bend_info(el)
    is_bend = false; rho = 0;
    if isfield(el,'h') && abs(el.h) > 1e-10
        is_bend = true;  rho = 1/abs(el.h);
    elseif isfield(el,'angle') && el.angle ~= 0 && el.length ~= 0
        is_bend = true;  rho = abs(el.length/el.angle);
    elseif isfield(el,'type')
        if any(strcmpi(el.type, {'dipole','sbend','rbend','bend','sector','rectangular'}))
            is_bend = true;
            if isfield(el,'h') && el.h ~= 0
                rho = 1/abs(el.h);
            else
                rho = 1.0;
            end
        end
    end
end


%% =====================================================================
%   Reversed-walk sum over drift segments between bends.
%   Mirrors the segment-treatment in ../../CalcCSRkick.m so that R56_seg
%   = R_cum(5,6) at the segment's start is the downstream-side R56.
%% =====================================================================
function [S_R56, S_one, n_seg] = local_drift_sums(beamline, sige, sigzf, sigz_in)
%   "drift" = any non-bend section between two bends (drift + quad + ...).
%   rho     = upstream bend's rho (original direction) = the bend we just hit
%             in the reversed walk.  Matches CalcCSRkick.m's convention.
    rev      = rev_beamline(beamline);
    R_cum    = eye(6);
    S_R56    = 0;
    S_one    = 0;
    n_seg    = 0;

    in_segment       = false;
    segment_start_R  = eye(6);
    segment_length   = 0;

    for k = 1:numel(rev)
        el = rev(k);
        [is_bend, rho] = local_bend_info(el);

        if is_bend
            if in_segment && segment_length > 0 && rho > 0
                R56_seg = segment_start_R(5,6);
                d_kern  = 1 - (sige/sigzf) * R56_seg;
                if d_kern > 0
                    C_seg = 1 / d_kern;
                else
                    C_seg = NaN;
                end
                phi_0 = (6 * sigz_in / abs(rho))^(1/3);
                phi_i = phi_0 / C_seg^(1/3);
                term  = 2 * segment_length / (abs(rho) * phi_i) + 1;
                if term <= 0
                    term = 1e-10;
                end
                lo    = log(term);
                S_one = S_one + C_seg * lo;
                S_R56 = S_R56 + R56_seg * C_seg * lo;
                n_seg = n_seg + 1;
            end
            in_segment     = false;
            segment_length = 0;
        else
            if ~in_segment
                in_segment      = true;
                segment_start_R = R_cum;
                segment_length  = el.length;
            else
                segment_length  = segment_length + el.length;
            end
        end

        R_cum = getElementMatrix(el) * R_cum;
    end
end


%% =====================================================================
%   Build lambda(z) on a uniform grid; normalized so that int lambda du = 1
%   with u = z / sigma_z.
%% =====================================================================
function [zg, lambda_grid] = local_lambda(z, n_grid, smooth_win, sigma_z)
    zmin   = min(z);
    zmax   = max(z);
    if zmax <= zmin
        error('local_lambda: degenerate z range.');
    end
    edges  = linspace(zmin, zmax, n_grid + 1).';
    dz     = edges(2) - edges(1);
    zg     = 0.5 * (edges(1:end-1) + edges(2:end));
    counts = histcounts(z, edges).';

    if smooth_win > 1
        smooth_counts = smoothdata(counts, 'gaussian', smooth_win);
        smooth_counts = smooth_counts * sum(counts) / sum(smooth_counts);
    else
        smooth_counts = counts;
    end

    N_tot       = sum(counts);
    lambda_grid = smooth_counts / N_tot * sigma_z / dz;
end
