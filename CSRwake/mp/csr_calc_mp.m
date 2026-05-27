function out = csr_calc_mp(z, delta, Q_total, beamline, gamma_rel, T566, varargin)
% CSR_CALC_MP  Multi-particle steady-state CSR + longitudinal transport,
%              driven by a lattice (struct array) rather than R56-line file.
%
%   out = csr_calc_mp(z, delta, Q_total, beamline, gamma_rel, T566, ...
%                     NAME, VALUE, ...)
%
%   The function walks through `beamline`, identifies bends from element.h
%   (or element.type), and accumulates the 6x6 transfer matrix. R56(s) is
%   taken as R_cumulative(5,6) at every substep inside each bend. Bend
%   radii rho are read directly from element.h.
%
%   Required
%     z         : [N x 1] particle z [m]
%     delta     : [N x 1] particle relative energy deviation
%     Q_total   : total bunch charge [C]
%     beamline  : struct array of lattice elements (compatible with
%                 getElementMatrix from ../../code). Each element needs at
%                 least .type and .length; bends additionally .h (=1/rho)
%                 or .angle.
%     gamma_rel : Lorentz factor
%     T566      : second-order longitudinal dispersion of the line [m].
%                 Pass [] to auto-compute it from the lattice via
%                 nonlinear_calculator(beamline, 5, 6, 6).
%
%   Optional name-value pairs
%     'sige'        ([])     override sige; default = std(delta)
%     'sigzf'       ([])     override sigzf; default = sigz0 + sige*r56_tot
%     'smooth_win'  (5)      Gaussian smoothing window for binned I(z)
%     'plot_flag'   (1)      diagnostic plots
%     'n_bends'     (0)      if >0, only use first n_bends bends
%     'dR56'        (0)      added to R56 in the transport step only
%     'n_grid'      (200)    grid points for lambda(z) build
%     'interp'      ('linear') interpolation method for Gamma -> particles
%     'n_substeps'  (200)    substeps per bend for the line integrals
%
%   Output: table with one row per particle, columns
%     z_out_noCSR, z_out_wCSR, delta_noCSR, delta_wCSR.
%   out.Properties.UserData carries Q_total, sige_used, R56_file, r56_tot,
%   sigma_z0, sigma_zf, S_R56, S_one, n_grid, smooth_win, plus diagnostic
%   grid arrays z_grid, I_grid, lambda_grid, Gamma_grid.

    % ---------- options ----------
    p = inputParser;
    p.addParameter('sige',       []);
    p.addParameter('sigzf',      []);
    p.addParameter('smooth_win', 5);
    p.addParameter('plot_flag',  1);
    p.addParameter('n_bends',    0);
    p.addParameter('dR56',       0);
    p.addParameter('n_grid',     200);
    p.addParameter('interp',     'linear');
    p.addParameter('n_substeps', 500);    % R56(s) sample points per dipole
                                          % (matches CalcCSRkick.m's bend integration)
    p.addParameter('drift_csr',  true);   % include drift / transient CSR (inter-bend)
    p.parse(varargin{:});
    opt = p.Results;

    % ---------- 0) sanitize inputs ----------
    if ~isstruct(beamline) || isempty(beamline)
        error('csr_calc_mp: beamline must be a non-empty struct array.');
    end

    % Auto-compute T566 from the FORWARD beamline if user passed [].
    if isempty(T566)
        E0_eV = gamma_rel * 0.511e6;
        T566  = nonlinear_calculator(beamline, 5, 6, 6, 1e-3, E0_eV);
        fprintf('csr_calc_mp: T566 auto-computed = %.6g m\n', T566);
    end

    beamline=rev_beamline(beamline);
    z     = z(:);
    delta = delta(:);
    if numel(z) ~= numel(delta)
        error('csr_calc_mp: z and delta must have the same length.');
    end
    if numel(z) < 10
        error('csr_calc_mp: too few particles (%d).', numel(z));
    end
    n_bad_z = nnz(~isfinite(z));
    n_bad_d = nnz(~isfinite(delta));
    if n_bad_z > 0 || n_bad_d > 0
        error(['csr_calc_mp: inputs contain non-finite values ' ...
               '(%d in z, %d in delta). Filter NaN/Inf upstream.'], ...
              n_bad_z, n_bad_d);
    end
    N_part = numel(z);

    % ---------- 1) sigz0, sige (particle-level, equal weights) ----------
    sigz0 = std(z);

    if isempty(opt.sige) || (isscalar(opt.sige) && isnan(opt.sige))
        sige = std(delta);
    else
        sige = opt.sige;
    end

    % ---------- 2) walk lattice for R56_total ----------
    R56_file = local_walk_total_R56(beamline);
    r56_tot  = R56_file + opt.dR56;

    % ---------- 3) rescale sige if dR56 ~= 0 (same as old version) -------
    if R56_file ~= 0 && opt.dR56 ~= 0
        sige = sige * (R56_file + opt.dR56) / R56_file;
    end

    % ---------- 4) sigzf (auto or override) -----------------------------
    %   Real post-line bunch length: include both R56*delta and T566*delta^2
    %   in the no-CSR transport, then std() over the full distribution.
    z_transported = z + r56_tot * delta + T566 * delta.^2;
    sigzf_auto    = std(z_transported);

    if isempty(opt.sigzf) || (isscalar(opt.sigzf) && isnan(opt.sigzf))
        sigzf = sigzf_auto;
        if sigzf <= 0
            error(['csr_calc_mp: sigma_zf <= 0 (got %.3g m). Past full ' ...
                   'compression -- check signs of sige, R56_file, dR56, T566.'], sigzf);
        end
    else
        sigzf = max(opt.sigzf, sigzf_auto);
    end

    % ---------- 5) line integrals over the bends via lattice walk -------
    [S_R56, S_one, n_bends_used] = local_lattice_integrals( ...
        beamline, sige, sigzf, opt.n_bends, opt.n_substeps);

    if n_bends_used == 0
        warning('csr_calc_mp: no bends found in beamline; S_R56 = S_one = 0.');
    end

    % ---------- 6) lambda(z) from particles + Gamma on grid -------------
    [z_grid, I_grid, lambda_grid, Gamma_grid, dDelta_grid, sigma_z_used] = ...
        local_csr_wake_mp(z, Q_total, gamma_rel, sigz0, sigzf, ...
                          opt.n_grid, opt.smooth_win);

    % ---------- 7) interpolate Gamma kick onto each particle ------------
    dDelta_part = interp1(z_grid, dDelta_grid, z, opt.interp, 0);

    % ---------- 8a) per-particle ss-CSR kicks ---------------------------
    ddelta_ss = dDelta_part * S_one;
    dz_ss     = dDelta_part * S_R56;

    % ---------- 8b) per-particle dr-CSR kicks (transient, drifts) ------
    if opt.drift_csr
        [S_R56_dr, S_one_dr, n_seg_dr] = ...
            local_drift_sums(beamline, sige, sigzf, sigz0);

        r_e_   = 2.8179403262e-15;
        qe_    = 1.602176634e-19;
        N_elec = Q_total / qe_;
        k_csrD = 2 * r_e_ * N_elec / (gamma_rel * sigzf);

        lambda_part = interp1(z_grid, lambda_grid, z, opt.interp, 0);

        ddelta_dr = k_csrD * lambda_part * S_one_dr;
        dz_dr     = k_csrD * lambda_part * S_R56_dr;
    else
        S_R56_dr = 0; S_one_dr = 0; n_seg_dr = 0;  k_csrD = 0;
        ddelta_dr = zeros(N_part, 1);
        dz_dr     = zeros(N_part, 1);
    end

    % ---------- 8c) combined per-particle CSR kicks --------------------
    ddelta_csr = ddelta_ss - ddelta_dr;
    dz_csr     = dz_ss     - dz_dr;

    % ---------- 9) transport (per particle) -----------------------------
    z_out_noCSR = z + r56_tot * delta + T566 * delta.^2;
    z_out_wCSR  = z_out_noCSR + dz_csr;
    delta_noCSR = delta;
    delta_wCSR  = delta + ddelta_csr;

    % ---------- 10) output table ----------------------------------------
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
        'sigma_z_u',       sigma_z_used, ...
        'S_R56',           S_R56,         ...   % ss-CSR scalar
        'S_one',           S_one,         ...   % ss-CSR scalar
        'n_bends_used',    n_bends_used,  ...
        'drift_csr',       opt.drift_csr, ...
        'S_R56_dr',        S_R56_dr,      ...   % dr-CSR scalar
        'S_one_dr',        S_one_dr,      ...   % dr-CSR scalar
        'n_segments_used', n_seg_dr,      ...
        'k_csrD',          k_csrD,        ...   % dr-CSR prefactor
        'n_grid',          opt.n_grid,    ...
        'smooth_win',      opt.smooth_win,...
        'z_grid',          z_grid,        ...
        'I_grid',          I_grid,        ...
        'lambda_grid',     lambda_grid,   ...
        'Gamma_grid',      Gamma_grid );

    % ---------- 11) plot ------------------------------------------------
    if opt.plot_flag
        figure('Name','csr\_calc\_mp','Color','w','Position',[100 100 1100 700]);

        subplot(2,2,1);
        plot(z_grid*1e6, I_grid, 'k-', 'LineWidth',1.4); grid on;
        xlabel('z  [\mum]'); ylabel('I(z)  [A]');
        title(sprintf('Binned + smoothed current (n\\_grid=%d, win=%d)', ...
                       opt.n_grid, opt.smooth_win));

        subplot(2,2,2);
        plot(z_grid/sigma_z_used, Gamma_grid, 'm-', 'LineWidth',1.4); grid on;
        xlabel('u = z/\sigma_z'); ylabel('\Gamma(u)');
        title('CSR kernel \Gamma(u)');

        subplot(2,2,3);
        plot(z_out_noCSR*1e6, delta_noCSR*1e3, '.', ...
             'Color',[.55 .75 1.0], 'MarkerSize',2); hold on;
        plot(z_out_wCSR *1e6, delta_wCSR *1e3, '.', ...
             'Color',[1.0 .45 .45], 'MarkerSize',2);
        grid on;
        xlabel('z_{out}  [\mum]'); ylabel('\delta_{out}  [10^{-3}]');
        legend('no CSR','with CSR','Location','best');
        title(sprintf('End-of-line phase space  (R_{56}=%.3g m, T_{566}=%.3g m)', ...
                       r56_tot, T566));

        subplot(2,2,4);
        plot(z*1e6, dz_csr*1e6,     'k.', 'MarkerSize',2); hold on;
        plot(z*1e6, ddelta_csr*1e3, 'r.', 'MarkerSize',2);
        grid on;
        xlabel('z_0  [\mum]');
        legend('dz_{csr}  [\mum]','d\delta_{csr}  [10^{-3}]', 'Location','best');
        title('Per-particle CSR kicks');
    end
end


%% =====================================================================
%%   Local helper: R56_total from a forward walk through the lattice
%% =====================================================================
function R56_total = local_walk_total_R56(beamline)
    R_cum = eye(6);
    for k = 1:numel(beamline)
        R_cum = getElementMatrix(beamline(k)) * R_cum;
    end
    R56_total = R_cum(5,6);
end


%% =====================================================================
%%   Local helper: classify a beamline element as a bend; return |rho|
%% =====================================================================
function [is_bend, rho] = local_bend_info(el)
    is_bend = false;
    rho     = 0;
    if isfield(el, 'h') && abs(el.h) > 1e-10
        is_bend = true;
        rho = 1 / abs(el.h);
    elseif isfield(el, 'angle') && el.angle ~= 0 && el.length ~= 0
        is_bend = true;
        rho = abs(el.length / el.angle);
    elseif isfield(el, 'type')
        if any(strcmpi(el.type, {'dipole','sbend','rbend','bend','sector','rectangular'}))
            is_bend = true;
            if isfield(el, 'h') && el.h ~= 0
                rho = 1 / abs(el.h);
            else
                rho = 1.0;
            end
        end
    end
end


%% =====================================================================
%%   Local helper: S_R56 and S_one by walking the lattice
%%   Inside each bend, substep with the element matrix and accumulate
%%   kernel = 1 / (1 - (sige/sigzf)*R56_acc(s))^{4/3} weighted by 1/|rho|^{2/3}
%% =====================================================================
function [S_R56, S_one, n_bends_used] = local_lattice_integrals( ...
        beamline, sige, sigzf, n_bends_limit, n_substeps)

    R_cum = eye(6);
    S_one = 0;
    S_R56 = 0;
    n_bends_used = 0;

    for k = 1:numel(beamline)
        el = beamline(k);
        [is_bend, rho] = local_bend_info(el);

        if is_bend && rho > 0 && el.length > 0
            % Skip past first n_bends_limit if specified
            if n_bends_limit > 0 && n_bends_used >= n_bends_limit
                R_cum = getElementMatrix(el) * R_cum;
                continue;
            end
            n_bends_used = n_bends_used + 1;

            ns = max(n_substeps, 2);
            ds = el.length / ns;

            el_step = el; el_step.length = ds;
            M_step  = getElementMatrix(el_step);
            el_half = el; el_half.length = ds/2;
            M_half  = getElementMatrix(el_half);

            inv_rho23 = 1 / rho^(2/3);

            R_cur = R_cum;
            for i = 1:ns
                R_mid    = M_half * R_cur;
                R56_at_s = R_mid(5,6);
                d        = 1 - (sige/sigzf) * R56_at_s;
                if d <= 0
                    kern = NaN;
                else
                    kern = 1 / d^(4/3);
                end
                contrib = kern * inv_rho23 * ds;
                S_one   = S_one + contrib;
                S_R56   = S_R56 + R56_at_s * contrib;
                R_cur   = M_step * R_cur;
            end
            R_cum = R_cur;
        else
            R_cum = getElementMatrix(el) * R_cum;
        end
    end
end


%% =====================================================================
%%   Local helper: build lambda(z) from particles, compute Gamma on grid
%% =====================================================================
function [zg, I_smooth, lambda, Gamma, dDelta_rho23, sigma_z] = ...
    local_csr_wake_mp(z, Q_total, gamma_rel, sigz0, sigzf, n_grid, smooth_win)

    re = 2.8179403262e-15;
    c  = 2.99792458e8;
    qe = 1.602176634e-19;

    N_part = numel(z);
    zmin   = min(z);
    zmax   = max(z);
    if zmax <= zmin
        error('local_csr_wake_mp: degenerate z range.');
    end

    edges = linspace(zmin, zmax, n_grid + 1).';
    dz    = edges(2) - edges(1);
    zg    = 0.5 * (edges(1:end-1) + edges(2:end));

    counts = histcounts(z, edges).';

    q_per_macro = Q_total / N_part;
    I_raw       = q_per_macro .* counts .* c ./ dz;

    if smooth_win > 1
        I_smooth = smoothdata(I_raw, 'gaussian', smooth_win);
        scale    = trapz(zg, I_raw) / trapz(zg, I_smooth);
        I_smooth = I_smooth * scale;
    else
        I_smooth = I_raw;
    end

    sigma_z =1.1 * sigz0;

    W      = trapz(zg, I_smooth);
    u      = zg / sigma_z;
    Qc     = W;
    lambda = I_smooth * sigma_z / Qc;

    dlam_du = gradient(lambda, u);
    Nu      = numel(u);
    Nt      = 400;
    Gamma   = zeros(Nu, 1);
    for k = 2:Nu
        Tk   = (u(k) - u(1))^(2/3);
        tt   = linspace(0, Tk, Nt).';
        uhat = u(k) - tt.^(3/2);
        fp   = interp1(u, dlam_du, uhat, 'linear', 0);
        Gamma(k) = -3 * (3/2) * trapz(tt, fp);
    end

    Q_T          = Qc / c;
    N_elec       = Q_T / qe;
    prefac       = 2 * re * N_elec / (3^(4/3) * gamma_rel * sigzf^(4/3));
    dDelta_rho23 = prefac * Gamma;
end


%% =====================================================================
%%   Local helper: dr-CSR sums (drift / inter-bend contribution)
%%   Walks the (already-reversed) beamline.  The "drift" is any section
%%   between two bends (drifts + quads + ...), and the rho used is the
%%   UPSTREAM bend's rho in original-line order = the bend just hit in
%%   the reversed walk.  Matches the convention of CalcCSRkick.m.
%% =====================================================================
function [S_R56, S_one, n_seg] = local_drift_sums(beamline_rev, sige, sigzf, sigz_in)
    R_cum = eye(6);
    S_R56 = 0;
    S_one = 0;
    n_seg = 0;

    in_segment      = false;
    segment_start_R = eye(6);
    segment_length  = 0;

    for k = 1:numel(beamline_rev)
        el = beamline_rev(k);
        [is_bend, rho] = local_bend_info(el);

        if is_bend
            % Process the just-accumulated segment using THIS bend's rho
            % (upstream bend in original-line order).
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
            % Any non-bend element (drift, quad, sext, marker, ...).
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
