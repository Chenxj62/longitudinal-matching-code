%% k(sigma_t) compare: Eq.(2) (scaled Hankel) vs Python resistive-wall approx
clear; close all; clc;

%% ---------------- Pipe / material parameters ----------------
b     = 25e-3;        % pipe radius [m] (50 mm diameter -> 25 mm radius)
sigma = 5.8e7;        % conductivity [S/m] (Cu)

eps0  = 8.854187817e-12;
mu0   = 4*pi*1e-7;
c     = 299792458;

%% ---------------- sigma_t scan (fs) ----------------
% 1000 fs -> 20 fs
Nsig = 80;
sigma_fs = logspace(log10(20), log10(1000), Nsig);
sigma_fs = fliplr(sigma_fs);          % 1000 -> 20
sigma_t  = sigma_fs * 1e-15;          % [s]

%% ---------------- omega grid for integration ----------------
% Gaussian weighting exp(-(sigma_t*omega)^2) kills tail above ~ few/sigma_t.
% Use omega_max based on smallest sigma_t (20 fs).
sigma_t_min = min(sigma_t);
omega_max = max(1e15, 80/sigma_t_min);      % [rad/s]  (80/sigma is conservative)
omega_min = 1;                               % avoid 0 for logspace; we add 0 separately
Nw = 20000;

omega = [0; logspace(log10(omega_min), log10(omega_max), Nw).'];  % column

%% ---------------- Re{Z||(omega)} from Eq.(2) ----------------
ReZ_eq2 = calcReZ_eq2_scaled(omega, b, sigma, eps0, mu0, c);  % [Ohm/m]
ReZ_eq2(~isfinite(ReZ_eq2)) = NaN;

%% ---------------- Re{Z||(omega)} from Python approx ----------------
% Python: R_wall(f) = (1/(2*pi*b)) * sqrt(mu0*(2*pi*f)/(2*sigma))
% Convert to omega: omega = 2*pi*f => R_wall(omega) = (1/(2*pi*b)) * sqrt(mu0*omega/(2*sigma))
ReZ_py = (1./(2*pi*b)) .* sqrt(mu0 .* omega ./ (2*sigma));    % [Ohm/m]
ReZ_py(omega==0) = 0;

%% ---------------- integrate k(sigma_t) for both ----------------
k_eq2 = zeros(size(sigma_t));
k_py  = zeros(size(sigma_t));

for ii = 1:numel(sigma_t)
    st = sigma_t(ii);
    W  = exp(-(st*omega).^2);

    k_eq2(ii) = (1/pi) * trapz(omega, ReZ_eq2 .* W);
    k_py(ii)  = (1/pi) * trapz(omega, ReZ_py  .* W);
end

%% ---------------- plot k comparison ----------------
figure('Position',[100 100 400 320]);
loglog(3e11*sigma_fs*1e-15, abs(k_eq2)*5e-11*5e-11*4.5e6, 'LineWidth', 1.6); hold on;
loglog(3e11*sigma_fs*1e-15, abs(k_py)*5e-11*5e-11*4.5e6,  '--', 'LineWidth', 1.6);
grid on; box on;
set(gca,'XDir','reverse');  % show 1000 -> 20 fs
xlabel('\sigma_z  (mm)');
ylabel('Power loss[W/m]');
legend('Analysis', '\lambda<<1', 'Location','best');
%title('k(\sigma_t) comparison: Eq.(2) vs Python resistive-wall approximation');

% %% ---------------- optional: inspect ReZ(omega) comparison ----------------
% % Plot ReZ vs frequency (Hz) over the same omega range (log-log).
% f = omega/(2*pi);
% figure('Position',[120 120 900 520]);
% semilogx(f(2:end), abs(ReZ_eq2(2:end)), 'LineWidth', 1.2); hold on;   % skip 0
% semilogx(f(2:end), abs(ReZ_py(2:end)),  '--', 'LineWidth', 1.2);
% grid on; box on;
% xlabel('f = \omega/2\pi (Hz)');
% ylabel('|Re\{Z_{||}\}| (Ohm/m)');
% legend('Eq.(2)', 'Python approx', 'Location','best');
% title('Re\{Z_{||}\} magnitude comparison (log-log)');

%% print a few values
fprintf("omega_max used = %.3e rad/s\n", omega_max);
fprintf("sigma_t(fs)        k_eq2            k_py           ratio(k_py/k_eq2)\n");
for jj = [1, round(Nsig/2), Nsig]
    fprintf("%10.3f   % .6e   % .6e   % .6e\n", sigma_fs(jj), k_eq2(jj), k_py(jj), k_py(jj)/k_eq2(jj));
end

%% ================= Local function: Eq.(2) stable Re{Z} =================
function ReZ = calcReZ_eq2_scaled(omega, b, sigma, eps0, mu0, c)
    omega = omega(:);
    ReZ = nan(size(omega));

    % omega=0 limit: ReZ -> 0 (resistive-wall ~ sqrt(omega))
    ReZ(omega==0) = 0;

    idx = omega > 0;
    w   = omega(idx);

    delta  = sqrt(2 ./ (sigma*mu0.*w));   % skin depth [m]
    lambda = (1+1i) ./ delta;             % [1/m]
    x1     = lambda .* b;

    % scaled Hankel to avoid overflow
    H0s = besselh(0, 1, x1, 1);
    H1s = besselh(1, 1, x1, 1);

    ratio = H1s ./ H0s;

    term2 = (w./(lambda*c) + (lambda*c)./w) .* ratio - (b.*w)./(2*c);
    Z2    = (-1i) ./ (2*pi*eps0*b*c .* term2);   % [Ohm/m]

    ReZ(idx) = real(Z2);
end
