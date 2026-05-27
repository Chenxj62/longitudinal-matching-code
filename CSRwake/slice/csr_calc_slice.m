function out = csr_calc_slice(z_slice, delta_slice, I_slice, beamline, ...
                              gamma_rel, T566, varargin)
% CSR_CALC_SLICE  Slice-based steady-state CSR + longitudinal transport,
%                 driven by a lattice (struct array).
%
%   out = csr_calc_slice(z_slice, delta_slice, I_slice, beamline, ...
%                        gamma_rel, T566, NAME, VALUE, ...)
%
%   Slice-level analog of csr_calc_mp.m. Identical physical model:
%     - walks `beamline`, identifies bends from element.h / element.type
%     - accumulates the 6x6 transfer matrix R_cum
%     - inside each bend, samples R56(s) = R_cum(5,6) at n_substeps midpoints
%     - line integrals S_R56, S_one with kernel 1/(1 - (sige/sigzf)*R56(s))^{4/3}
%     - Gamma(u) HFP integral built from the slice I(z); per-slice CSR kicks.
%
%   The difference from csr_calc_mp.m is granularity: lambda(z) is built
%   directly from the input slice arrays (no histogramming), and the CSR
%   kicks are evaluated at each slice center (no interpolation onto particles).
%   Output has one row per input slice.
%
%   Required
%     z_slice     : [Ns x 1] slice center z [m]
%     delta_slice : [Ns x 1] slice-mean delta
%     I_slice     : [Ns x 1] slice current [A]
%     beamline    : struct array of lattice elements (compatible with
%                   getElementMatrix). Bends need .h or .angle.
%     gamma_rel   : Lorentz factor
%     T566        : second-order longitudinal dispersion [m]. Pass [] to
%                   auto-compute via nonlinear_calculator(beamline,5,6,6).
%
%   Optional name-value pairs
%     'sige'        ([])     override sige; default = current-weighted std of delta_slice
%     'sigzf'       ([])     override sigzf; default = sigz0 + sige*r56_tot
%     'smooth_win'  (5)      Gaussian smoothing window for I(z)
%     'plot_flag'   (1)      diagnostic plots
%     'n_bends'     (0)      if >0, only use first n_bends bends
%     'dR56'        (0)      added to R56 in the transport step only
%     'n_substeps'  (500)    R56(s) sample points per dipole
%
%   The beamline is always walked in reverse (CalcCSRkick convention), so
%   R_cum(5,6) at any substep is the residual R56 from that point to the
%   end of the original line.
%
%   Output: table with one row per slice, columns
%     z0, z_out_noCSR, z_out_wCSR, delta_noCSR, delta_wCSR, dIdz_smooth
%
%   out.Properties.UserData carries the resolved scalars:
%     sige_used, R56_file, r56_tot, sigma_z0, sigma_zf, S_R56, S_one,
%     n_bends_used, smooth_win, plus arrays z_grid, I_grid, lambda_grid,
%     Gamma_grid for diagnostics.

    % ---------- options ----------
    p = inputParser;
    p.addParameter('sige',        []);
    p.addParameter('sigzf',       []);
    p.addParameter('smooth_win',  5);
    p.addParameter('plot_flag',   1);
    p.addParameter('n_bends',     0);
    p.addParameter('dR56',        0);
    p.addParameter('n_substeps',  500);
    p.addParameter('drift_csr',   true);    % include drift / transient CSR
    p.parse(varargin{:});
    opt = p.Results;

    % ---------- 0) sanitize inputs ----------
    z_slice     = z_slice(:);
    delta_slice = delta_slice(:);
    I_slice     = I_slice(:);

    if numel(z_slice) ~= numel(delta_slice) || numel(z_slice) ~= numel(I_slice)
        error('csr_calc_slice: z_slice, delta_slice, I_slice must have same length.');
    end
    if numel(z_slice) < 5
        error('csr_calc_slice: too few slices (%d).', numel(z_slice));
    end
    if any(~isfinite(z_slice)) || any(~isfinite(delta_slice)) || any(~isfinite(I_slice))
        error('csr_calc_slice: inputs contain non-finite values. Filter upstream.');
    end
    if ~isstruct(beamline) || isempty(beamline)
        error('csr_calc_slice: beamline must be a non-empty struct array.');
    end

    [z_slice, ix] = sort(z_slice);
    delta_slice   = delta_slice(ix);
    I_slice       = I_slice(ix);

    % Auto-compute T566 from the FORWARD beamline if user passed [].
    if isempty(T566)
        E0_eV = gamma_rel * 0.511e6;
        T566  = nonlinear_calculator(beamline, 5, 6, 6, 1e-3, E0_eV);
        fprintf('csr_calc_slice: T566 auto-computed = %.6g m\n', T566);
    end

    beamline = rev_beamline(beamline);


    % ---------- 1) sigz0, sige (current-weighted moments of slice data) -
    [zbar0, sigz0] = local_weighted_moments(z_slice,     I_slice);
    if isempty(opt.sige) || (isscalar(opt.sige) && isnan(opt.sige))
        [~, sige] = local_weighted_moments(delta_slice, I_slice);
    else
        sige = opt.sige;
    end

    % ---------- 2) walk lattice for R56_total ---------------------------
    R56_file = local_walk_total_R56(beamline);
    r56_tot  = R56_file + opt.dR56;

    % ---------- 3) rescale sige if dR56 ~= 0 ----------------------------
    if R56_file ~= 0 && opt.dR56 ~= 0
        sige = sige * (R56_file + opt.dR56) / R56_file;
    end

    % ---------- 4) sigzf (auto or override) -----------------------------
    %   Real post-line bunch length: include both R56*delta and T566*delta^2
    %   in the no-CSR transport, then current-weighted std over slices.
    z_transported   = z_slice + r56_tot * delta_slice + T566 * delta_slice.^2;
    [~, sigzf_auto] = local_weighted_moments(z_transported, I_slice);

    if isempty(opt.sigzf) || (isscalar(opt.sigzf) && isnan(opt.sigzf))
        sigzf = sigzf_auto;
        if sigzf <= 0
            error(['csr_calc_slice: sigma_zf <= 0 (got %.3g m). Past full ' ...
                   'compression -- check signs of sige, R56_file, dR56, T566.'], sigzf);
        end
    else
        sigzf = max(opt.sigzf, sigzf_auto);
    end

    % ---------- 5) line integrals over the bends via lattice walk -------
    [S_R56, S_one, n_bends_used] = local_lattice_integrals( ...
        beamline, sige, sigzf, opt.n_bends, opt.n_substeps);

    if n_bends_used == 0
        warning('csr_calc_slice: no bends found in beamline; S_R56 = S_one = 0.');
    end

    % ---------- 6) lambda(z) on slice grid + Gamma ----------------------
    [I_smooth, lambda_grid, Gamma_grid, dDelta_grid, sigma_z_used] = ...
        local_csr_wake_slice(z_slice, I_slice, gamma_rel, sigz0, sigzf, ...
                             opt.smooth_win);

    % ---------- 7a) per-slice ss-CSR kicks ------------------------------
    ddelta_ss = dDelta_grid * S_one;
    dz_ss     = dDelta_grid * S_R56;

    % ---------- 7b) per-slice dr-CSR kicks (drift / inter-bend) --------
    if opt.drift_csr
        [S_R56_dr, S_one_dr, n_seg_dr] = ...
            local_drift_sums(beamline, sige, sigzf, sigz0);

        r_e_   = 2.8179403262e-15;
        qe_    = 1.602176634e-19;
        c_     = 2.99792458e8;
        Qc     = trapz(z_slice, I_smooth);    % = Q_total * c
        N_elec = (Qc / c_) / qe_;
        k_csrD = 2 * r_e_ * N_elec / (gamma_rel * sigzf);

        ddelta_dr = k_csrD * lambda_grid * S_one_dr;
        dz_dr     = k_csrD * lambda_grid * S_R56_dr;
    else
        S_R56_dr = 0; S_one_dr = 0; n_seg_dr = 0; k_csrD = 0;
        ddelta_dr = zeros(size(z_slice));
        dz_dr     = zeros(size(z_slice));
    end

    % ---------- 7c) combined per-slice CSR kicks -----------------------
    ddelta_csr = ddelta_ss - ddelta_dr;
    dz_csr     = dz_ss     - dz_dr;

    % ---------- 8) transport (per slice) --------------------------------
    z_out_noCSR = z_slice + r56_tot * delta_slice + T566 * delta_slice.^2;
    z_out_wCSR  = z_out_noCSR + dz_csr;
    delta_noCSR = delta_slice;
    delta_wCSR  = delta_slice + ddelta_csr;

    % ---------- 9) output table -----------------------------------------
    dIdz_smooth = gradient(I_smooth, z_slice);
    out = table(z_slice, z_out_noCSR, z_out_wCSR, ...
                delta_noCSR, delta_wCSR, dIdz_smooth, ...
                'VariableNames', {'z0','z_out_noCSR','z_out_wCSR', ...
                                  'delta_noCSR','delta_wCSR','dIdz_smooth'});
    out.Properties.UserData = struct( ...
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
        'smooth_win',      opt.smooth_win,...
        'z_grid',          z_slice,       ...
        'I_grid',          I_smooth,      ...
        'lambda_grid',     lambda_grid,   ...
        'Gamma_grid',      Gamma_grid );

    % ---------- 10) plot ------------------------------------------------
    if opt.plot_flag
        % --- Figure 1: phase space w/wo CSR + dI_in_smooth/ds ----------
        figure('Name','csr\_calc\_slice','Color','w','Position',[100 100 1100 450]);

        subplot(1,2,1);
        plot(z_out_noCSR*1e6, delta_noCSR*1e3, 'b.-', ...
             'LineWidth',1.0, 'MarkerSize',6); hold on;
        plot(z_out_wCSR *1e6, delta_wCSR *1e3, 'r.-', ...
             'LineWidth',1.0, 'MarkerSize',6);
        grid on;
        xlabel('z_{out}  [\mum]'); ylabel('\delta_{out}  [10^{-3}]');
        legend('no CSR','with CSR','Location','best');
        title(sprintf('End-of-line phase space  (R_{56}=%.3g m, T_{566}=%.3g m)', ...
                       r56_tot, T566));

        subplot(1,2,2);
        plot(z_slice*1e6, dIdz_smooth, 'm-', 'LineWidth',1.4); grid on;
        xlabel('z  [\mum]'); ylabel('dI_{smooth}/dz  [A/m]');
        title('Derivative of smoothed input current');

        % --- Figure 2: output current via slice Jacobian ----------------
        zp_noCSR  = gradient(z_out_noCSR, z_slice);
        zp_wCSR   = gradient(z_out_wCSR,  z_slice);
        I_out_noCSR = I_slice ./ max(abs(zp_noCSR), eps);
        I_out_wCSR  = I_slice ./ max(abs(zp_wCSR),  eps);

        figure('Name','csr\_calc\_slice: I_{out}','Color','w', ...
               'Position',[120 120 700 450]);
        plot(z_out_noCSR*1e6, I_out_noCSR, 'b-', 'LineWidth',1.4); hold on;
        plot(z_out_wCSR *1e6, I_out_wCSR,  'r-', 'LineWidth',1.4);
        grid on;
        xlabel('z_{out}  [\mum]'); ylabel('I_{out}  [A]');
        legend('no CSR','with CSR','Location','best');
        title(sprintf(['Output current (slice Jacobian)  ' ...
                       'peak: no CSR=%.1f A, w/ CSR=%.1f A'], ...
                       max(I_out_noCSR), max(I_out_wCSR)));
    end
end


%% =====================================================================
%%   Local helper: current-weighted (xbar, sigma) of a slice array
%% =====================================================================
function [xbar, sig] = local_weighted_moments(x, w)
    W = sum(w);
    if W <= 0
        xbar = mean(x); sig = std(x);
        return;
    end
    xbar = sum(x .* w) / W;
    sig  = sqrt(sum((x - xbar).^2 .* w) / W);
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
%%   Local helper: bend / non-bend classification, returns |rho|
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
%% =====================================================================
function [S_R56, S_one, n_bends_used] = local_lattice_integrals( ...
        beamline, sige, sigzf, n_bends_limit, n_substeps)

    R_cum        = eye(6);
    S_one        = 0;
    S_R56        = 0;
    n_bends_used = 0;

    for k = 1:numel(beamline)
        el = beamline(k);
        [is_bend, rho] = local_bend_info(el);

        if is_bend && rho > 0 && el.length > 0
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
%%   Local helper: slice-level CSR wake (delta'_csr * rho^{2/3})
%%   I_slice already lives on z_slice (no histogramming).
%% =====================================================================
function [Isl, lambda, Gamma, dDelta_rho23, sigma_z] = ...
    local_csr_wake_slice(zc, Isl_raw, gamma_rel, sigz0, sigzf, smooth_win)

    re = 2.8179403262e-15;
    c  = 2.99792458e8;
    qe = 1.602176634e-19;

    zc      = zc(:);
    Isl_raw = Isl_raw(:);

    if smooth_win > 1
        Isl = smoothdata(Isl_raw, 'gaussian', smooth_win);
        % preserve total charge
        scale = trapz(zc, Isl_raw) / trapz(zc, Isl);
        Isl   = Isl * scale;
    else
        Isl = Isl_raw;
    end

    sigma_z = 1.0 * sigz0;

    W      = trapz(zc, Isl);
    u      = zc / sigma_z;
    Qc     = W;
    lambda = Isl * sigma_z / Qc;          % int lambda du = 1

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
    dDelta_rho23 = prefac * Gamma;        % [m^-1/3]
end


%% =====================================================================
%%   Local helper: dr-CSR sums (drift / inter-bend contribution)
%%   Same algorithm as csr_calc_mp.m: walks the (already-reversed)
%%   beamline; for each inter-bend segment, uses the just-hit bend's rho
%%   in the reversed walk = upstream bend's rho in original direction.
%%   Matches the convention of CalcCSRkick.m.
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
