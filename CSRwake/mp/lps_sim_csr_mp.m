function out = lps_sim_csr_mp(p)
% LPS_SIM_CSR_MP  Multi-particle LPS pipeline.
%
%   out = lps_sim_csr_mp(p)
%
%   Stages:
%     0. read_dist(bunchfile)                     -> z0, delta0  (z centered)
%     1. Linac 1 (no-op on delta; sets E_ref1)
%     2. Arc 1 via csr_calc_mp                    -> z2, delta2
%     3. Linac 2 (per-particle RF chirp)          -> z_l2, delta_l2
%     4. Arc 2 via csr_calc_mp                    -> z_out, delta_out
%
%   Required fields of p
%     bunchfile        distribution file (ASTRA-style; cols 3=z, 6=p, 10=stat)
%     Q_total          total bunch charge [C]
%     E_ref0           initial reference energy [MeV]
%     A                cavity voltage shared by both linacs [MeV]
%     phi1, phi2       linac phases [deg]
%     f                RF frequency [Hz]
%     T566_2           arc-2 second-order R56 [m]
%     r56_arc2_file    R56(s) line file for arc 2
%     rho_list_arc2    bend radii for arc 2 [m]
%
%   Optional fields (defaults shown)
%     T566_1         (0)         arc-1 second-order R56 [m]
%     r56_arc1_file  ('r56')     R56(s) file for arc 1
%     rho_list_arc1  (1)         bend radii for arc 1 [m]
%     n_grid         (400)       grid points for lambda(z) build
%     smooth_win     (40)        Gaussian smoothing window for I(z)
%     n_bends_arc1   (0)         restrict arc-1 csr_calc_mp to first n bends
%     n_bends_arc2   (0)         restrict arc-2 csr_calc_mp to first n bends
%     dR56_1         (0)         dR56 added to arc-1 transport step only
%     dR56_2         (0)         dR56 added to arc-2 transport step only
%     gamma_rel      ([])        Lorentz factor for arc 2; auto if empty
%     sigzf_arc1     ([])        sigma_zf override for arc 1
%     sigzf_arc2     ([])        sigma_zf override for arc 2
%     s_col          (2)         s column in r56 files
%     R56_col        (3)         R56 column in r56 files
%     plot           (0)         1 to draw a diagnostic figure
%     verbose        (1)         1 to print stage diagnostics
%
%   Returned struct fields
%     z0, delta0                  initial particle arrays (z centered)
%     z2, delta2                  after arc 1 (with CSR)
%     z_l2, delta_l2              after linac 2
%     z_out, delta_out            final (with CSR)
%     E_ref1, E_ref2              reference energies
%     gamma_arc1, gamma_arc2      gamma used in each csr_calc_mp
%     sige_arc1, sige_arc2        sige used in each arc
%     R56_arc1, R56_arc2          R56_file used in each arc
%     N_part, Q_total
%     is_valid                    false if any stage produced non-finite particles
%     fail_stage                  'arc1' / 'arc2' / '' if ok

    % ---------- defaults ----------
    if ~isfield(p,'T566_1'),         p.T566_1         = 0;     end
    if ~isfield(p,'r56_arc1_file'),  p.r56_arc1_file  = 'r56'; end
    if ~isfield(p,'rho_list_arc1'),  p.rho_list_arc1  = 1;     end
    if ~isfield(p,'n_grid'),         p.n_grid         = 400;   end
    if ~isfield(p,'smooth_win'),     p.smooth_win     = 40;    end
    if ~isfield(p,'n_bends_arc1'),   p.n_bends_arc1   = 0;     end
    if ~isfield(p,'n_bends_arc2'),   p.n_bends_arc2   = 0;     end
    if ~isfield(p,'dR56_1'),         p.dR56_1         = 0;     end
    if ~isfield(p,'dR56_2'),         p.dR56_2         = 0;     end
    if ~isfield(p,'gamma_rel'),      p.gamma_rel      = [];    end
    if ~isfield(p,'sigzf_arc1'),     p.sigzf_arc1     = [];    end
    if ~isfield(p,'sigzf_arc2'),     p.sigzf_arc2     = [];    end
    if ~isfield(p,'s_col'),          p.s_col          = 2;     end
    if ~isfield(p,'R56_col'),        p.R56_col        = 3;     end
    if ~isfield(p,'plot'),           p.plot           = 0;     end
    if ~isfield(p,'verbose'),        p.verbose        = 1;     end
    if ~isfield(p,'L_cav_2'),        p.L_cav_2        = 0;     end

    c        = 2.99792458e8;
    k_rf     = 2*pi * p.f / c;
    phi1_rad = deg2rad(p.phi1);
    phi2_rad = deg2rad(p.phi2);

    %% Stage 0: load particles (z centered)
    [z0, delta0] = read_dist(p.bunchfile);
    N_part = numel(z0);

    %% Stage 1: Linac 1
    %   (delta unchanged in current setup; uncomment the cavity_long line
    %    when you want the Linac-1 RF chirp active)
    % [z1, delta1, E_ref1] = cavity_long(z0, delta0, p.E_ref0, p.A, p.phi1, p.f);
    E_ref1 = p.E_ref0 + p.A * cos(phi1_rad);
    delta1 = delta0;
    z1     = z0;

    %% Stage 2: Arc 1
    gamma_arc1 = E_ref1 / 0.510998950;
    args1 = build_csr_args(p, 1);

    arc1 = csr_calc_mp(z1, delta1, p.Q_total, p.r56_arc1_file, gamma_arc1, ...
                       p.rho_list_arc1, p.T566_1, args1{:});
    z2     = arc1.z_out_wCSR;
    delta2 = arc1.delta_wCSR;

    ud1 = arc1.Properties.UserData;
    n_bad_arc1 = nnz(~isfinite(z2)) + nnz(~isfinite(delta2));

    if n_bad_arc1 > 0
        out = make_invalid_out(p, N_part);
        out.E_ref1 = E_ref1; out.gamma_arc1 = gamma_arc1;
        out.z0 = z0; out.delta0 = delta0;
        out.fail_stage = 'arc1';
        return;
    end

    %% Stage 3: Linac 2 (thin-cavity full-cosine; same as lps_sim_csr.m)
    [z_l2, delta_l2, E_ref2] = cavity_long(z2, delta2, E_ref1, ...
                                           p.A, p.phi2, p.f);

    %% Stage 4: Arc 2
    if isempty(p.gamma_rel)
        gamma_arc2 = E_ref2 / 0.510998950;
    else
        gamma_arc2 = p.gamma_rel;
    end
    args2 = build_csr_args(p, 2);

    arc2 = csr_calc_mp(z_l2, delta_l2, p.Q_total, p.r56_arc2_file, gamma_arc2, ...
                       p.rho_list_arc2, p.T566_2, args2{:});
    z_out     = arc2.z_out_wCSR;
    delta_out = arc2.delta_wCSR;

    ud2 = arc2.Properties.UserData;
    n_bad_arc2 = nnz(~isfinite(z_out)) + nnz(~isfinite(delta_out));

    %% pack
    out.z0          = z0;
    out.delta0      = delta0;
    out.z2          = z2;
    out.delta2      = delta2;
    out.z_l2        = z_l2;
    out.delta_l2    = delta_l2;
    out.z_out       = z_out;
    out.delta_out   = delta_out;
    out.E_ref1      = E_ref1;
    out.E_ref2      = E_ref2;
    out.gamma_arc1  = gamma_arc1;
    out.gamma_arc2  = gamma_arc2;
    out.sige_arc1   = ud1.sige_used;
    out.sige_arc2   = ud2.sige_used;
    out.R56_arc1    = ud1.R56_file;
    out.R56_arc2    = ud2.R56_file;
    out.N_part      = N_part;
    out.Q_total     = p.Q_total;
    out.is_valid    = (n_bad_arc1 == 0) && (n_bad_arc2 == 0);
    out.fail_stage  = '';

    %% optional diagnostic plot (currents via histogram at every stage)
    if p.plot
        n_bins = max(200, round(p.n_grid));
        qpm    = p.Q_total / N_part;

        [zc_in,  I_in ] = histcur(z0,    n_bins, qpm, c);
        [zc_a1,  I_a1 ] = histcur(z2,    n_bins, qpm, c);
        [zc_l2,  I_l2 ] = histcur(z_l2,  n_bins, qpm, c);
        [zc_o,   I_o  ] = histcur(z_out, n_bins, qpm, c);

        figure('Color','w','Name','lps\_sim\_csr\_mp', ...
               'Position',[80 80 1100 800]);

        subplot(2,2,1);
        zg1 = ud1.z_grid;  lam1 = ud1.lambda_grid;
        zg2 = ud2.z_grid;  lam2 = ud2.lambda_grid;
        dlam1 = gradient(lam1, zg1);
        dlam2 = gradient(lam2, zg2);
        plot(zg1*1e6, dlam1, 'b-', 'LineWidth',1.4); hold on;
        plot(zg2*1e6, dlam2, 'r-', 'LineWidth',1.4); grid on;
        xlabel('z  [\mum]'); ylabel('d\lambda/dz  [1/m]');
        legend('Arc 1 in','Arc 2 in','Location','best');
        title('\lambda''(z) entering each arc CSR kernel');

        subplot(2,2,2);
        plot(z_out*1e6, delta_out*1e3, '.', ...
             'Color',[1.0 .45 .45], 'MarkerSize',2);
        grid on;
        xlabel('z_{out}  [\mum]'); ylabel('\delta_{out}  [10^{-3}]');
        title('Final phase space');

        subplot(2,2,3);
        plot(zc_o*1e6, I_o, 'r-', 'LineWidth',1.4); grid on;
        xlabel('z_{out}  [\mum]'); ylabel('I_{out}  [A]');
        title(sprintf('Final I_{out}   peak = %.1f A', max(I_o)));

        subplot(2,2,4);
        dd_arc2 = arc2.delta_wCSR - arc2.delta_noCSR;
        plot(z_l2*1e6, dd_arc2*1e3, 'r.', 'MarkerSize',2); grid on;
        xlabel('z (Arc-2 in)  [\mum]'); ylabel('\delta_{csr, arc2}  [10^{-3}]');
        title('Arc-2 CSR \delta kick');
    end
end


%% =====================================================================
function args = build_csr_args(p, arc_idx)
    if arc_idx == 1
        args = {'n_grid',     p.n_grid, ...
                'smooth_win', p.smooth_win, ...
                'n_bends',    p.n_bends_arc1, ...
                'dR56',       p.dR56_1, ...
                's_col',      p.s_col, ...
                'R56_col',    p.R56_col, ...
                'plot_flag',  0};
        if ~isempty(p.sigzf_arc1)
            args = [args, {'sigzf', p.sigzf_arc1}];
        end
    else
        args = {'n_grid',     p.n_grid, ...
                'smooth_win', p.smooth_win, ...
                'n_bends',    p.n_bends_arc2, ...
                'dR56',       p.dR56_2, ...
                's_col',      p.s_col, ...
                'R56_col',    p.R56_col, ...
                'plot_flag',  0};
        if ~isempty(p.sigzf_arc2)
            args = [args, {'sigzf', p.sigzf_arc2}];
        end
    end
end


%% =====================================================================
%   ASTRA-style distribution reader.
%     - Row 1 is the reference particle (absolute coords).
%     - Rows 2..N store z and p as deltas relative to row 1.
%   Unwind row 1, drop it, subtract mean(z) so the bunch sits around z=0,
%   then form delta = (p - <p>) / <p>.
%% =====================================================================
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

    keep = isfinite(z_abs) & isfinite(p_abs);
    z = z_abs(keep);
    p = p_abs(keep);

    z = z - mean(z);                 % center longitudinal coordinate

    pbar = mean(p);
    if ~isfinite(pbar) || pbar == 0
        error('read_dist: mean(p) = %g; check that column 6 is momentum.', pbar);
    end
    delta = (p - pbar) / pbar;
end


%% =====================================================================
function [zc, I] = histcur(z, nb, qpm, c)
    zmin  = min(z);
    zmax  = max(z);
    edges = linspace(zmin, zmax, nb + 1);
    dz    = edges(2) - edges(1);
    zc    = 0.5 * (edges(1:end-1) + edges(2:end));
    cnt   = histcounts(z, edges);
    I     = qpm * cnt * c / dz;
    zc    = zc(:);
    I     = I(:);
end


%% =====================================================================
function out = make_invalid_out(p, N_part)
    out = struct( ...
        'z0',         NaN, 'delta0',     NaN, ...
        'z2',         NaN, 'delta2',     NaN, ...
        'z_l2',       NaN, 'delta_l2',   NaN, ...
        'z_out',      NaN, 'delta_out',  NaN, ...
        'E_ref1',     NaN, 'E_ref2',     NaN, ...
        'gamma_arc1', NaN, 'gamma_arc2', NaN, ...
        'sige_arc1',  NaN, 'sige_arc2',  NaN, ...
        'R56_arc1',   NaN, 'R56_arc2',   NaN, ...
        'N_part',     N_part, ...
        'Q_total',    p.Q_total, ...
        'is_valid',   false, ...
        'fail_stage', '');
end
