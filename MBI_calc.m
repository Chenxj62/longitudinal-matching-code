function result = MBI_calc(beamline, beam, opts)
% MBI_calc  Microbunching instability gain — CSR + optics
%
% Self-contained Volterra integral solver. No globals, no MBIcode needed.
% Based on C.-Y. Tsai, PRAB 19, 114401 (2016).
%
% Usage:
%   result = MBI_calc(beamline, beam, opts)
%
% Inputs:
%   beamline  - struct array from createElement
%   beam      - struct:
%       .energy_MeV, .I_b, .emit_norm_x, .emit_norm_y  [MeV, A, um, um]
%       .betax0, .betay0, .alphax0, .alphay0             [m, m, -, -]
%       .sigma_delta, .chirp (optional, m^-1), .sigma_z  [-, m^-1, m]
%   opts      - struct (all optional):
%       .lambda_range     [um_start, um_end]  (default: auto)
%       .scan_num         (default: 50)
%       .mesh_num         (default: 400)
%       .iCSR_ss          steady-state CSR    (default: 1)
%       .iCSR_tr          transient CSR       (default: 0)
%       .iCSR_drift       CSR drift effect    (default: 0)
%       .issCSRpp         CSR shielding       (default: 0)
%       .full_pipe_height pipe height [cm]    (default: 1e50)
%       .itransLD         transverse LD       (default: 1)
%       .calc_energy_gain                     (default: 1)
%       .plot                                 (default: 0)

%% Defaults
if nargin < 3, opts = struct(); end
def = struct('lambda_range',[], 'scan_num',50, 'mesh_num',400, ...
    'iCSR_ss',1, 'iCSR_tr',0, 'iCSR_drift',0, ...
    'issCSRpp',0, 'full_pipe_height',1e50, ...
    'itransLD',1, 'calc_energy_gain',1, 'plot',0);
for f = fieldnames(def)'
    if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
end
if ~isfield(beam,'chirp'), beam.chirp = 0; end

%% Constants
me_MeV   = 0.511;
re_cgs   = 2.81794e-13;     % classical electron radius [cm]
I_A      = 17045;            % Alfven current [A]
A_csr    = 3^(-1/3) * gamma(2/3) * (1i*sqrt(3) - 1);

gamma0   = beam.energy_MeV / me_MeV;
ex_cm    = beam.emit_norm_x * 1e-4 / gamma0;   % geometric emittance [cm]
ey_cm    = beam.emit_norm_y * 1e-4 / gamma0;
bx0      = beam.betax0 * 100;    % [cm]
by0      = beam.betay0 * 100;
ax0      = beam.alphax0;
ay0      = beam.alphay0;
sd       = beam.sigma_delta;
chirp_cm = beam.chirp / 100;     % [cm^-1]
nb       = beam.I_b / (re_cgs * I_A);

% ======================================================================
%  Step 1: Walk beamline — R(0→s), gamma(s), dipole table
% ======================================================================
nelem = length(beamline);
s_si  = zeros(nelem+1, 1);
R_all = cell(nelem+1, 1);
g_arr = zeros(nelem+1, 1);

R_all{1} = eye(6);  g_arr(1) = gamma0;
R_cum = eye(6);  cur_g = gamma0;

dip_ent = []; dip_ext = []; dip_rho = [];

for ie = 1:nelem
    [Re, go] = getElementMatrix(beamline(ie), cur_g);
    R_cum = Re * R_cum;  cur_g = go;
    s_si(ie+1)   = s_si(ie) + beamline(ie).length;
    R_all{ie+1}  = R_cum;
    g_arr(ie+1)  = cur_g;
    if strcmpi(beamline(ie).type,'dipole') && abs(beamline(ie).h) > 1e-10
        dip_ent(end+1) = ie;  dip_ext(end+1) = ie+1;  %#ok
        dip_rho(end+1) = 1/beamline(ie).h;             %#ok
    end
end
Npts = nelem + 1;

% ======================================================================
%  Step 2: Extract R5x, convert SI → CGS
% ======================================================================
s_cm = s_si * 100;
R51=zeros(Npts,1); R52=R51; R53=R51; R54=R51; R55=R51; R56=R51;
for ip = 1:Npts
    R = R_all{ip};
    R51(ip)=R(5,1); R52(ip)=R(5,2); R53(ip)=R(5,3); R54(ip)=R(5,4);
    R55(ip)=R(5,5); R56(ip)=R(5,6);
end
% CGS scaling
R52 = R52*100;  R54 = R54*100;  R56 = R56*100;
% HK sign convention
R51 = -R51;  R52 = -R52;  R53 = -R53;  R54 = -R54;  R56 = -R56;
% Velocity bunching R56
R56 = R56 + (s_cm - s_cm(1)) ./ g_arr.^2;
% Compression
C = abs(1 ./ (R55 - chirp_cm * R56));

% ======================================================================
%  Step 3: Dipole table [cm]
% ======================================================================
N_D = length(dip_ent);
dip = struct('N', N_D, 's1', zeros(N_D,1), 's2', zeros(N_D,1), ...
    'rho', zeros(N_D,1), 'rho_x', zeros(N_D,1), 'rho_y', zeros(N_D,1));
for id = 1:N_D
    dip.s1(id)    = s_cm(dip_ent(id));
    dip.s2(id)    = s_cm(dip_ext(id));
    dip.rho(id)   = dip_rho(id) * 100;        % [cm]
    dip.rho_x(id) = dip_rho(id) * 100;
    dip.rho_y(id) = 1e50;
end

% ======================================================================
%  Step 4: Mesh and precompute
% ======================================================================
N  = opts.mesh_num;
s0 = s_cm(1);  s1 = s_cm(end) - 0.1;
sm = linspace(s0, s1, N)';
h  = (s1 - s0) / N;

Cm   = interp1(s_cm, C,   sm, 'linear','extrap');
R51m = interp1(s_cm, R51, sm, 'linear','extrap');
R52m = interp1(s_cm, R52, sm, 'linear','extrap');
R53m = interp1(s_cm, R53, sm, 'linear','extrap');
R54m = interp1(s_cm, R54, sm, 'linear','extrap');
R55m = interp1(s_cm, R55, sm, 'linear','extrap');
R56m = interp1(s_cm, R56, sm, 'linear','extrap');

% Bending radius at each mesh point
rho_m = ones(N,1) * 1e50;
for id = 1:N_D
    mask = (sm >= dip.s1(id)) & (sm <= dip.s2(id));
    rho_m(mask) = dip.rho(id);
end

% ======================================================================
%  Step 5: Auto wavelength range
% ======================================================================
if isempty(opts.lambda_range)
    R56end = abs(R56m(end));
    if R56end < 1e-4
        lopt = 3 * mean(abs(R56m)) * sd * 1e4;
    else
        lopt = 3 * R56end * sd * 1e4;
    end
    lam_start = max(0.2*lopt, 0.01);
    lam_end   = max(50*lopt, 100);
else
    lam_start = opts.lambda_range(1);
    lam_end   = opts.lambda_range(2);
end
Nscan = opts.scan_num;
iEcalc = opts.calc_energy_gain;

% ======================================================================
%  Step 6: Volterra solver — wavelength scan
% ======================================================================
lam_arr = linspace(lam_start, lam_end, Nscan);
Gf_dd = zeros(1,Nscan);  G_vs_s = zeros(N,Nscan);
if iEcalc
    Gf_de = zeros(1,Nscan); Gf_ed = zeros(1,Nscan);
    Gf_ee = zeros(1,Nscan);
end

coeff0 = 1i * re_cgs * nb / gamma0;

fprintf('MBI_calc: %d wavelengths [%.2f, %.2f] um ...\n', Nscan, lam_start, lam_end);

for iscan = 1:Nscan
    lam_cm = lam_arr(iscan) * 1e-4;
    k = 2*pi / lam_cm;

    % --- Landau damping at mesh (for G0 vectors) ---
    LDx = 0.5*k^2*ex_cm/bx0 * Cm.^2 .* (bx0^2*(R51m-R52m*ax0/bx0).^2 + R52m.^2);
    LDy = 0.5*k^2*ey_cm/by0 * Cm.^2 .* (by0^2*(R53m-R54m*ay0/by0).^2 + R54m.^2);
    LDd = 0.5*k^2*sd^2 * Cm.^2 .* R56m.^2;
    if opts.itransLD, LD = LDx+LDy+LDd; else, LD = LDd; end

    % --- G0 vectors ---
    n_1k0 = 1.0;  e_1k0 = 0.0;
    G0_k  = n_1k0 * exp(-LD);
    G0_kp = -1i*k * Cm .* R56m * e_1k0 .* exp(-LD);
    E0_k  = -1i*k * Cm .* R56m * sd^2 * n_1k0 .* exp(-LD);
    E0_kp = (1 - k^2*Cm.^2.*R56m.^2*sd^2) * e_1k0 .* exp(-LD);

    % --- Build kernel matrix ---
    K_mat  = zeros(N);
    ML_mat = zeros(N);

    for j = 1:(N-1)
        sj = sm(j);
        rj = rho_m(j);
        Cj = Cm(j);

        % === CSR impedance Z(k*C(s'), s') ===
        Z = compute_Z_csr(k, Cj, sj, rj, dip, gamma0, opts, A_csr);
        if Z == 0, continue; end

        ii = (j+1):N;
        fac = 1; if j == 1, fac = 0.5; end

        % Modified transport R5x_mod(s,s') = C(s)*R5x(s) - C(s')*R5x(s')
        dR51 = Cm(ii).*R51m(ii) - Cj*R51m(j);
        dR52 = Cm(ii).*R52m(ii) - Cj*R52m(j);
        dR53 = Cm(ii).*R53m(ii) - Cj*R53m(j);
        dR54 = Cm(ii).*R54m(ii) - Cj*R54m(j);
        dR56 = Cm(ii).*R56m(ii) - Cj*R56m(j);

        % R56 transport (symplectic)
        R56t = R56m(ii) - R56m(j)*R55m(ii)/R55m(j) ...
             + R51m(j)*R52m(ii) - R51m(ii)*R52m(j) ...
             + R53m(j)*R54m(ii) - R53m(ii)*R54m(j);

        % Landau damping in kernel
        LDx_k = 0.5*k^2*ex_cm/bx0 * (bx0^2*(dR51-dR52*ax0/bx0).^2 + dR52.^2);
        LDy_k = 0.5*k^2*ey_cm/by0 * (by0^2*(dR53-dR54*ay0/by0).^2 + dR54.^2);
        LDd_k = 0.5*k^2*sd^2 * dR56.^2;
        if opts.itransLD, LD_k = LDx_k+LDy_k+LDd_k; else, LD_k = LDd_k; end

        % kernel_mod (density-density)
        K_mat(ii, j) = fac * coeff0*k .* Cm(ii) .* Cj .* Z .* R56t .* exp(-LD_k);

        if iEcalc
            % kernel_mod_M (energy coupling M)
            K_M = fac * re_cgs*nb/gamma0 * k^2 .* Cm(ii).^2 .* Cj ...
                * sd^2 .* Z .* dR56 .* R56t .* exp(-LD_k);
            % kernel_mod_L (energy coupling L)
            K_L = fac * re_cgs*nb/gamma0 .* Cj .* Z .* exp(-LD_k);
            ML_mat(ii, j) = K_M - K_L;
        end
    end

    % --- Solve Volterra ---
    I_K = eye(N) - h * K_mat;
    gkd = I_K \ G0_k;
    G_c_dd = gkd / max(abs(G0_k(1)), 1e-30);
    Gf_dd(iscan) = abs(G_c_dd(end));
    G_vs_s(:, iscan) = abs(G_c_dd);

    if iEcalc
        gkp = I_K \ G0_kp;
        ekd = h * ML_mat * (I_K \ G0_k)  + E0_k;
        ekp = h * ML_mat * (I_K \ G0_kp) + E0_kp;
        Gf_ed(iscan) = abs(gkp(end) / max(abs(e_1k0),1e-30));
        Gf_de(iscan) = abs(ekd(end) / max(abs(G0_k(1)),1e-30));
        Gf_ee(iscan) = abs(ekp(end) / max(abs(e_1k0),1e-30));
    end

    if mod(iscan, max(floor(Nscan/10),1))==0 || iscan==Nscan
        fprintf('  [%d/%d] lam=%.2f um, Gf=%.4f\n', iscan, Nscan, lam_arr(iscan), Gf_dd(iscan));
    end
end

% ======================================================================
%  Step 7: Output
% ======================================================================
result.lambda   = lam_arr;
result.s        = sm / 100;       % [m]
result.Gf_dd    = Gf_dd;
result.G_vs_s   = G_vs_s;
if iEcalc
    result.Gf_de    = Gf_de;
    result.Gf_ed    = Gf_ed;
    result.Gf_ee    = Gf_ee;
    result.Gf_total = abs(Gf_dd + Gf_ed);
    result.Ef_total = abs(Gf_de + Gf_ee);
end

% ======================================================================
%  Step 8: Plot
% ======================================================================
if opts.plot && Nscan >= 2
    figure; plot(lam_arr, Gf_dd, 'b-', 'LineWidth', 2);
    xlabel('\lambda (\mum)'); ylabel('Density gain'); grid on;
    title('MBI: density-to-density');
    if iEcalc
        figure; plot(lam_arr, result.Gf_total, 'k-', 'LineWidth', 2);
        xlabel('\lambda (\mum)'); ylabel('Total density gain'); grid on;
        title('MBI: total density gain');
    end
    figure; surf(lam_arr, result.s, G_vs_s);
    xlabel('\lambda (\mum)'); ylabel('s (m)'); zlabel('|G|');
    title('MBI gain map'); shading interp; colorbar;
end

[pk, pidx] = max(Gf_dd);
fprintf('MBI_calc done. Peak gain = %.4f at lambda = %.2f um\n', pk, lam_arr(pidx));
end


% ======================================================================
%  CSR impedance: all effects combined
% ======================================================================
function Z = compute_Z_csr(k, C_sp, s_p, rho_p, dip, egamma, opts, A_csr)
% Returns total CSR impedance Z(k*C(s'), s') in CGS
% k: wavenumber [cm^-1], C_sp: compression at s', s_p: position [cm]
% rho_p: bending radius at s' [cm], dip: dipole table struct

    kC = k * C_sp;
    Z = 0;
    in_dipole = abs(rho_p) < 1e20;

    % --- Steady-state CSR ---
    if opts.iCSR_ss && in_dipole
        if opts.issCSRpp
            % Shielding: check threshold wavenumber
            h_pipe = opts.full_pipe_height;
            k_th = 2*pi * sqrt(abs(rho_p) / h_pipe^3);
            if kC < k_th
                Z = Z + csr_shielding(kC, h_pipe, rho_p);
            else
                Z = Z + csr_ss(kC, rho_p, A_csr);
            end
        else
            Z = Z + csr_ss(kC, rho_p, A_csr);
        end
    end

    % --- Transient CSR ---
    if opts.iCSR_tr && in_dipole
        Z = Z + csr_tr(kC, s_p, rho_p, dip);
    end

    % --- CSR drift ---
    if opts.iCSR_drift && ~in_dipole
        Z = Z + csr_drift(kC, s_p, dip, egamma);
    end
end


% ======================================================================
%  Steady-state CSR (ultra-relativistic, free-space)
% ======================================================================
function Z = csr_ss(kC, rho, A_csr)
    Z = -1i * kC^(1/3) * A_csr * abs(rho)^(-2/3);
end


% ======================================================================
%  Transient CSR
%  Inside dipole, within formation length: entrance transient
%  Ref: Stupakov, PRST-AB 10 (2007); Tsai (2016) Eq.
% ======================================================================
function Z = csr_tr(kC, s_p, rho, dip)
    local_s = find_local_s_func(s_p, dip);
    if local_s == 0, Z = 0; return; end
    if local_s < 0.05, local_s = 0.05; end

    lam_tmp = 2*pi / kC;
    formation_length = (24 * lam_tmp * rho^2)^(1/3);

    if local_s < formation_length
        z_L = local_s^3 / (24 * rho^2);
        mu  = kC * z_L;
        t1  = -4/local_s * exp(-1i*4*mu);
        t2  = 4/(3*local_s) * (1i*mu)^(1/3) * upper_igamma(-1/3, 1i*mu);
        Z   = t1 + t2;
    else
        Z = 0;
    end
end


% ======================================================================
%  CSR drift (CER wake after dipole exit)
%  Ref: R. Bosch, PRST-AB 11, 090702 (2008)
% ======================================================================
function Z = csr_drift(kC, s_p, dip, egamma)
    if dip.N == 0, Z = 0; return; end

    % Find all upstream dipoles
    upstream_Lb  = [];
    upstream_rho = [];
    downstream_s = [];
    for m = 1:dip.N
        if s_p > dip.s2(m)
            upstream_Lb(end+1)  = abs(dip.s2(m) - dip.s1(m));  %#ok
            upstream_rho(end+1) = dip.rho(m);                  %#ok
            downstream_s(end+1) = s_p - dip.s2(m);             %#ok
        end
    end

    if isempty(upstream_Lb), Z = 0; return; end

    Lb  = upstream_Lb(:);
    rho = upstream_rho(:);
    ds  = downstream_s(:);

    % Case D: Bosch Eq.(24), PRST-AB 10, 050701 (2007)
    r = length(ds);
    tmp03 = zeros(r,1);
    for m = 1:r
        if ds(m) > abs(rho(m))^(2/3) * (2*pi/kC)^(1/3)
            d_eff = min(ds(m), egamma^2/kC);
            tmp03(m) = 2 / d_eff;
        end
    end

    % Case C
    tmp04 = exp(-1i*kC*Lb.^2./(6*rho.^2).*(Lb + 3*ds));
    tmp05 = Lb + 2*ds;
    tmp06 = -4*tmp04 ./ tmp05;

    Z = tmp03(end) + tmp06(end);
end


% ======================================================================
%  CSR shielding (parallel plates)
%  Ref: Agoh & Yokoya, PRSTAB 7, 054403 (2004)
% ======================================================================
function Z = csr_shielding(kC, h_pipe, rho)
    p_max = 200;
    rho_abs = abs(rho);
    coeff = 8*pi^2/h_pipe * (2/(kC*rho_abs))^(1/3);
    s = 0;
    for p = 0:p_max
        beta_p = (2*p+1)*pi/h_pipe * (rho_abs/(2*kC^2))^(1/3);
        s = s + F0_func(beta_p);
    end
    Z = coeff * s;
end

function f = F0_func(b)
    if b > 10
        t1 = b/(2*pi) * exp(-4/3*b^3);
        t2 = 1 + 1/(24*b^3) + 1/(1152*b^6);
        t3 = -3/(16*pi*b^5);
        t4 = 1 + 105/32/b^6;
        f  = t1*t2 + 1i*t3*t4;
    else
        b2 = b^2;
        t1 = airy(1,b2).*airy(1,b2) - 1i*airy(1,b2).*airy(3,b2);
        t2 = b2.*airy(0,b2).*airy(0,b2) - 1i*b2.*airy(0,b2).*airy(2,b2);
        f  = t1 + t2;
    end
end


% ======================================================================
%  Helper: distance from dipole entrance
% ======================================================================
function ds = find_local_s_func(s, dip)
    ds = 0;
    for m = 1:dip.N
        if s >= dip.s1(m) && s <= dip.s2(m)
            ds = s - dip.s1(m);
        end
    end
end


% ======================================================================
%  Upper incomplete gamma: Gamma(a,z) for complex z
%  Replaces mfun('GAMMA', a, z) without Symbolic Math Toolbox
%  Algorithm from Numerical Recipes (gammaincc/gammalnz in MBIcode)
% ======================================================================
function val = upper_igamma(a, z)
    % Gamma(a,z) = Gamma(a) * (1 - P(a,z))
    ga = exp(gammalnz_local(a));
    val = ga * (1 - gammaincc_local(z, a));
end

function b = gammaincc_local(x, a)
    if x == 0, b = 0; return; end
    lnga = gammalnz_local(a);
    if abs(x) < abs(a)+1
        % Series expansion
        ap = a; s = 1/a; del = s;
        for n = 1:500
            ap = ap + 1;
            del = x * del / ap;
            s = s + del;
            if abs(del) < 10*eps*abs(s), break; end
        end
        b = s * exp(-x + a*log(x) - lnga);
    else
        % Continued fraction
        a0 = 1; a1 = x; b0 = 0; b1 = 1;
        fac = 1; g = b1; gold = b0;
        n = 1;
        while abs(g - gold) >= 10*eps*abs(g)
            gold = g;
            ana = n - a;
            a0 = (a1 + a0*ana) * fac;
            b0 = (b1 + b0*ana) * fac;
            anf = n * fac;
            a1 = x * a0 + anf * a1;
            b1 = x * b0 + anf * b1;
            fac = 1 / a1;
            g = b1 * fac;
            n = n + 1;
            if n > 500, break; end
        end
        b = 1 - exp(-x + a*log(x) - lnga) * g;
    end
end

function f = gammalnz_local(z)
    cof = [76.18009172947146; -86.50532032941677; 24.01409824083091; ...
           -1.231739572450155; 0.1208650973866179e-2; -0.5395239384953e-5];
    zp = z + (1:6)';
    ser = sum(cof ./ zp) + 1.000000000190015;
    tmp = z + 5.5 - (z + 0.5) * log(z + 5.5);
    f = -tmp + log(2.5066282746310005 * ser / z);
end
