function results = BBU_dipole_tracking(p)
% BBU_DIPOLE_TRACKING  Time-domain dipole-mode BBU bunch tracking
%
% Translated from UniBBU/BunchTracking.py (dipole modes only).
%
% Usage:
%   results = BBU_dipole_tracking()       % default parameters
%   results = BBU_dipole_tracking(p)      % custom parameters via struct
%
% Parameters (struct fields, all optional — defaults in parentheses):
%   machine_index     : machine realization index (0)
%   active_modes      : 0-based mode indices to activate ([27, 28])
%   f0                : fundamental RF frequency [Hz] (1.3e9)
%   Nc                : number of cavities (4)
%   Na                : number of cavity passes (16)
%   NmTotal           : total modes per cavity in file (128)
%   Nb                : number of bunches to track (100000)
%   I_beam            : average beam current [A] (0.855)
%   initial_V_dip     : initial dipole HOM voltage seed [V] (1.0)
%   long_motion       : include longitudinal offsets in timing (false)
%   output_dir        : directory for output files ('results')
%   max_volt_lines    : max rows in Voltage.txt (100)
%   output_time_units : time unit for output files (1e-3 = ms)
%   machines_dir      : base path for input files ('machines')
%
% Returns struct with: decrement, volt_time, volt_amp, X_final,
%                      f_hom, active_modes, stable

%% ========== CONSTANTS ==========
Ee = 0.5109989461e6;   % electron rest energy [eV]
c  = 2.99792458e8;     % speed of light [m/s]

%% ========== PARAMETERS ==========
if nargin < 1, p = struct(); end

machine_index     = dval(p, 'machine_index', 0);
active_modes      = dval(p, 'active_modes', [27, 28]);
f0                = dval(p, 'f0', 1.3e9);
Nc                = dval(p, 'Nc', 4);
Na                = dval(p, 'Na', 16);
NmTotal           = dval(p, 'NmTotal', 128);
Nb                = dval(p, 'Nb', 100000);
I_beam            = dval(p, 'I_beam', 0.855);
initial_V_dip     = dval(p, 'initial_V_dip', 1.0);
long_motion       = dval(p, 'long_motion', false);
output_dir        = dval(p, 'output_dir', 'results');
max_volt_lines    = dval(p, 'max_volt_lines', 100);
output_time_units = dval(p, 'output_time_units', 1e-3);
machines_dir      = dval(p, 'machines_dir', 'machines');

%% ========== DERIVED QUANTITIES ==========
t0 = 1.0 / f0;
q  = I_beam * t0;      % bunch charge [C]
Nm = length(active_modes);

%% ========== READ TRANSPORT MATRICES ==========
% File format per arc (Na blocks):
%   line 1 : energy [eV]
%   line 2 : arc length [m]
%   lines 3-8 : 6x6 transport matrix (one row per line)

matrix_file = fullfile(machines_dir, 'matrices', sprintf('%04d.txt', machine_index));
fid = fopen(matrix_file, 'r');
if fid < 0
    error('Cannot open matrix file: %s', matrix_file);
end

Ea = zeros(Na, 1);     % energy per arc [eV]
La = zeros(Na, 1);     % arc length [m]
Ta = cell(Na, 1);      % 6x6 transport matrix per arc

for ia = 1:Na
    Ea(ia) = fscanf(fid, '%f', 1);
    La(ia) = fscanf(fid, '%f', 1);
    T = zeros(6, 6);
    for row = 1:6
        T(row, :) = fscanf(fid, '%f', 6)';
    end
    Ta{ia} = T;
end
fclose(fid);

Va   = c * sqrt(1.0 - (Ee ./ Ea).^2);   % velocity per arc [m/s]
Zcav = cumsum(La);                        % cumulative cavity positions [m]

%% ========== READ CAVITY HOM PARAMETERS ==========
% File format: Nc cavities x NmTotal lines each.
% Columns: index  order  freq[Hz]  Q  R/Q[Ohm/m^2]  theta[deg]  [df  dQQ]
% We keep only dipole modes (order == 1) whose index is in active_modes.

cavity_file = fullfile(machines_dir, 'cavities', sprintf('%04d.txt', machine_index));
fid = fopen(cavity_file, 'r');
if fid < 0
    error('Cannot open cavity file: %s', cavity_file);
end

f_hom   = zeros(Nc, Nm);   % HOM frequency [Hz]
Q_hom   = zeros(Nc, Nm);   % quality factor
RoQd    = zeros(Nc, Nm);   % R/Q [Ohm/m^2]
Theta_h = zeros(Nc, Nm);   % polarization angle [rad]

for ic = 1:Nc
    im_found = 0;
    for j = 1:NmTotal
        line = fgetl(fid);             % read one line (handles any column count)
        A = sscanf(line, '%f')';
        if ismember(A(1), active_modes) && A(2) == 1   % dipole mode check
            im_found = im_found + 1;
            f_hom(ic, im_found)   = A(3);
            Q_hom(ic, im_found)   = A(4);
            RoQd(ic, im_found)    = A(5);
            Theta_h(ic, im_found) = A(6) * pi / 180.0;
        end
    end
    if im_found ~= Nm
        warning('Cavity %d: found %d active dipole modes, expected %d', ...
            ic, im_found, Nm);
    end
end
fclose(fid);

w_hom   = 2.0 * pi * f_hom;           % angular frequency [rad/s]
Gamma_h = w_hom ./ (2.0 * Q_hom);     % half-bandwidth [1/s]

%% ========== INITIALIZE BUNCHES ==========
BunchSpacing    = Va(1) * t0;
X               = zeros(Nb, 6);            % [x x' y y' z delta]
Z               = -(0:Nb-1)' * BunchSpacing;
BunchUpdateTime = zeros(Nb, 1);

BunchesToCheck        = (Nb + 1) * ones(Na, 1);
BunchesToCheck(1)     = 1;

%% ========== INITIALIZE CAVITY VOLTAGES ==========
% Seed with real value (matches Python: InitialModeAmplitudes[1] for dipole).
% The oscillation exp((-Gamma+i*w)*dt) generates the imaginary part naturally.
V_cav            = initial_V_dip * ones(Nc, Nm);
CavityUpdateTime = zeros(Nc, 1);

%% ========== GROWTH RATE ACCUMULATORS ==========
sumt    = zeros(Nc, Nm);
sumt2   = zeros(Nc, Nm);
sumlnV  = zeros(Nc, Nm);
sumtlnV = zeros(Nc, Nm);
countN  = zeros(Nc, Nm);
decrement = zeros(Nc, Nm);

t_half = 0.50 * Nb * t0;
t_3qtr = 0.75 * Nb * t0;

%% ========== OUTPUT SETUP ==========
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
DataPointsPerWrite = max(1, floor(Nb * Na / Nc / max_volt_lines));
LinesCount = zeros(Nc, 1);

% Per-cavity voltage buffers (fixes interleaving bug in original)
max_rows  = ceil(Nb * Na / Nc / DataPointsPerWrite) + 10;
volt_time = cell(Nc, 1);
volt_amp  = cell(Nc, 1);
volt_row  = zeros(Nc, 1);
for ic_init = 1:Nc
    volt_time{ic_init} = zeros(max_rows, 1);
    volt_amp{ic_init}  = zeros(max_rows, Nm);
end

%% ========== MAIN TRACKING LOOP ==========
fprintf('Tracking %d bunches through %d arcs, I = %.4f A ...\n', Nb, Na, I_beam);
t       = 0.0;
t_start = tic;
total_turns = Na * Nb;
next_report = 0.25;

for turn = 1:total_turns

    % ---- Progress report ----
    if turn / total_turns >= next_report
        fprintf('  %.0f%% done (%.1f s elapsed)\n', next_report * 100, toc(t_start));
        next_report = next_report + 0.25;
    end

    % ---- Find the next bunch-cavity arrival event ----
    dtmin  = inf;
    ia_evt = 0;
    ib_evt = 0;
    for ia = 1:Na
        ib = BunchesToCheck(ia);
        if ib <= Nb
            dt = (Zcav(ia) - Z(ib)) / Va(ia) - (t - BunchUpdateTime(ib));
            if dt < dtmin
                dtmin  = dt;
                ia_evt = ia;
                ib_evt = ib;
            end
        end
    end

    % ---- Advance global time ----
    t  = t + dtmin;
    ia = ia_evt;
    ib = ib_evt;

    % ---- Transport bunch through arc ----
    X(ib, :) = (Ta{ia} * X(ib, :)')';
    Z(ib)    = Zcav(ia);
    BunchUpdateTime(ib) = t;

    BunchesToCheck(ia) = BunchesToCheck(ia) + 1;
    if ia < Na && ib == 1
        BunchesToCheck(ia + 1) = 1;
    end

    % ---- Map arc index to cavity index (periodic) ----
    ic = mod(ia - 1, Nc) + 1;

    % ---- Time since last cavity update ----
    if ~long_motion
        dt_cav = t - CavityUpdateTime(ic);
    else
        dt_cav = t + X(ib, 5) / c - CavityUpdateTime(ic);
    end

    % ---- Relativistic momentum*c [eV]:  pc = E*v/c ----
    if ia < Na
        pc = 0.5 * (Ea(ia) * Va(ia) + Ea(ia + 1) * Va(ia + 1)) / c;
    else
        pc = Ea(ia) * Va(ia) / c;
    end

    % ---- Dipole HOM interaction ----
    for im = 1:Nm
        costh = cos(Theta_h(ic, im));
        sinth = sin(Theta_h(ic, im));

        % 1) Free oscillation + beam excitation
        V_cav(ic, im) = V_cav(ic, im) ...
            * exp((-Gamma_h(ic, im) + 1i * w_hom(ic, im)) * dt_cav) ...
            + RoQd(ic, im) * q * c * (X(ib, 1) * costh + X(ib, 3) * sinth);

        % 2) Transverse kick from imaginary part of cavity voltage
        kick = imag(V_cav(ic, im)) / (2.0 * pc);
        X(ib, 2) = X(ib, 2) + kick * costh;
        X(ib, 4) = X(ib, 4) + kick * sinth;
    end

    % ---- Update cavity clock ----
    if ~long_motion
        CavityUpdateTime(ic) = t;
    else
        CavityUpdateTime(ic) = t + X(ib, 5) / c;
    end

    % ---- Growth rate: accumulate log|V| for linear regression ----
    if t > t_half && t < t_3qtr
        lnV = log(abs(V_cav(ic, :)));
        sumt(ic, :)    = sumt(ic, :)    + t;
        sumt2(ic, :)   = sumt2(ic, :)   + t^2;
        sumlnV(ic, :)  = sumlnV(ic, :)  + lnV;
        sumtlnV(ic, :) = sumtlnV(ic, :) + t * lnV;
        countN(ic, :)  = countN(ic, :)  + 1;
    end

    % ---- Store downsampled voltage (per cavity) ----
    if mod(LinesCount(ic), DataPointsPerWrite) == 0
        volt_row(ic) = volt_row(ic) + 1;
        volt_time{ic}(volt_row(ic))    = t / output_time_units;
        volt_amp{ic}(volt_row(ic), :)  = abs(V_cav(ic, :));
    end
    LinesCount(ic) = LinesCount(ic) + 1;

end  % main loop

fprintf('Tracking done. Elapsed: %.2f s\n', toc(t_start));

%% ========== COMPUTE GROWTH RATES ==========
fprintf('\n--- Growth rates (decrement) ---\n');
any_unstable = false;
for ic = 1:Nc
    for im = 1:Nm
        n = countN(ic, im);
        if n > 1
            decrement(ic, im) = ...
                (n * sumtlnV(ic, im) - sumlnV(ic, im) * sumt(ic, im)) / ...
                (n * sumt2(ic, im)   - sumt(ic, im)^2);
        end
        if decrement(ic, im) > 0
            status = 'UNSTABLE';
            any_unstable = true;
        else
            status = 'stable';
        end
        fprintf('  Cavity %d  Mode %-3d  f = %9.4f MHz  decrement = %+.4e 1/s  [%s]\n', ...
            ic, active_modes(im), f_hom(ic, im) / 1e6, decrement(ic, im), status);
    end
end
if any_unstable
    fprintf('\nResult: UNSTABLE\n');
else
    fprintf('\nResult: stable\n');
end

%% ========== TRIM OUTPUT ARRAYS ==========
for ic = 1:Nc
    volt_time{ic} = volt_time{ic}(1:volt_row(ic));
    volt_amp{ic}  = volt_amp{ic}(1:volt_row(ic), :);
end

%% ========== SAVE OUTPUT FILES ==========
% Voltage.txt — merged: [t_c1 |V_c1_m1| |V_c1_m2| ... t_c2 |V_c2_m1| ...]
min_rows = min(volt_row);
merged = zeros(min_rows, Nc * (1 + Nm));
for ic = 1:Nc
    col0 = (ic - 1) * (1 + Nm) + 1;
    merged(:, col0)             = volt_time{ic}(1:min_rows);
    merged(:, col0+1:col0+Nm)  = volt_amp{ic}(1:min_rows, :);
end
writematrix(merged, fullfile(output_dir, 'Voltage.txt'), 'Delimiter', '\t');

% Decrement.txt
writematrix(decrement, fullfile(output_dir, 'Decrement.txt'), 'Delimiter', '\t');

% FinalCoords.txt
writematrix(X, fullfile(output_dir, 'FinalCoords.txt'), 'Delimiter', '\t');

fprintf('\nResults saved to %s/\n', output_dir);

%% ========== PLOT VOLTAGE EVOLUTION ==========
figure('Name', 'BBU Dipole HOM Voltage Evolution');
for ic = 1:Nc
    for im = 1:Nm
        subplot(Nc, Nm, (ic - 1) * Nm + im);
        semilogy(volt_time{ic}, volt_amp{ic}(:, im), 'b-', 'LineWidth', 1);
        xlabel(sprintf('Time [%.0e s]', output_time_units));
        ylabel('|V| [V]');
        title(sprintf('Cav %d  Mode %d  (%.4f MHz)', ...
            ic, active_modes(im), f_hom(ic, im) / 1e6));
        grid on;
    end
end

%% ========== RETURN RESULTS ==========
results.decrement    = decrement;
results.volt_time    = volt_time;
results.volt_amp     = volt_amp;
results.X_final      = X;
results.f_hom        = f_hom;
results.active_modes = active_modes;
results.stable       = ~any_unstable;

end  % function BBU_dipole_tracking


%% ========== LOCAL HELPER ==========
function val = dval(s, name, default)
    if isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
