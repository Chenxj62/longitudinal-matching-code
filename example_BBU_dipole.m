% example_BBU_dipole.m
% Example: generate synthetic ERL machine data and run dipole BBU tracking
%
% This sets up a simple 2-cavity, 4-pass energy recovery linac and
% tracks 20000 bunches to determine BBU stability.
%
% Machine layout (2 cavities, 1 recirculation):
%   Source -> [Cav1] -> drift -> [Cav2] -> recirculation arc -> [Cav1] -> drift -> [Cav2] -> dump
%   pass 1: inject at 5 MeV, arrive at Cav1
%   pass 2: accelerated to 25 MeV, arrive at Cav2
%   pass 3: accelerated to 45 MeV, recirculate back to Cav1
%   pass 4: decelerated to 25 MeV, arrive at Cav2, then dump

clear; clc; close all;

%% ========== MACHINE PARAMETERS ==========
Nc = 2;        % number of cavities
Na = 4;        % number of cavity passes (arcs)
NmTotal = 4;   % total HOM entries per cavity in file

% Energy [eV] at each cavity pass
Ea = [5e6; 25e6; 45e6; 25e6];

% Arc length [m] between consecutive cavities
La = [0.5; 2.0; 15.0; 2.0];

%% ========== BUILD TRANSPORT MATRICES ==========
T = cell(Na, 1);

% Pass 1: short drift from source to cavity 1 (L = 0.5 m)
T{1} = eye(6);
T{1}(1,2) = La(1);   T{1}(3,4) = La(1);

% Pass 2: FODO cell, cavity 1 -> cavity 2 (mu = pi/4, L = 2 m)
mu = pi/4;  L = La(2);
T{2} = fodo_matrix(mu, L);

% Pass 3: recirculation arc, cavity 2 -> cavity 1 (mu = pi/3, L = 15 m)
mu = pi/3;  L = La(3);
T{3} = fodo_matrix(mu, L);

% Pass 4: FODO cell, cavity 1 -> cavity 2 (same as pass 2)
T{4} = T{2};

%% ========== CAVITY HOM DATA ==========
% 4 modes per cavity: 2 monopole (order=0, inactive) + 2 dipole (order=1, active)
% Columns: index  order  freq[Hz]  Q  R/Q[Ohm/m^2]  theta[deg]  df  dQQ
%
% Dipole mode parameters are based on TESLA 1.3 GHz cavity modes 27 & 28.
base_modes = [
    0  0  1.280000e+09  1.000000e+06  1.000000e-03    0.0   0  0   % monopole (not used)
    1  0  1.300000e+09  5.000000e+05  2.000000e-03    0.0   0  0   % monopole (not used)
    2  1  1.733559e+09  1.090000e+04  7.935390e+04  -85.9   0  0   % dipole: TESLA mode 27
    3  1  1.733838e+09  1.660000e+04  7.894105e+04    4.0   0  0   % dipole: TESLA mode 28
];

%% ========== WRITE INPUT FILES ==========
machines_dir = 'machines';
mat_dir = fullfile(machines_dir, 'matrices');
cav_dir = fullfile(machines_dir, 'cavities');
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end
if ~exist(cav_dir, 'dir'), mkdir(cav_dir); end

% --- Transport matrices ---
fid = fopen(fullfile(mat_dir, '0000.txt'), 'w');
for ia = 1:Na
    fprintf(fid, '%.6e\n', Ea(ia));
    fprintf(fid, '%.6e\n', La(ia));
    for row = 1:6
        fprintf(fid, '%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', T{ia}(row,:));
    end
end
fclose(fid);

% --- Cavity HOMs ---
fid = fopen(fullfile(cav_dir, '0000.txt'), 'w');
for ic = 1:Nc
    freq_shift = (ic - 1) * 1.0e6;   % 1 MHz spread between cavities
    for j = 1:NmTotal
        m = base_modes(j, :);
        fprintf(fid, '%d\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', ...
            m(1), m(2), m(3) + freq_shift, m(4), m(5), m(6), m(7), m(8));
    end
end
fclose(fid);

fprintf('Input files written to %s/\n', machines_dir);

%% ========== PRINT MACHINE SUMMARY ==========
fprintf('\n--- Machine summary ---\n');
fprintf('  Cavities: %d,  Passes: %d\n', Nc, Na);
fprintf('  Energies [MeV]:');
fprintf(' %.1f', Ea / 1e6);
fprintf('\n');
fprintf('  Arc lengths [m]:');
fprintf(' %.1f', La);
fprintf('\n');
fprintf('  Recirculation T12 = %.2f m  (arc 3)\n', T{3}(1,2));
fprintf('  Active dipole modes: 2, 3\n');
fprintf('  Mode 2: f = %.4f MHz, Q = %.0f, R/Q = %.0f Ohm/m^2\n', ...
    base_modes(3,3)/1e6, base_modes(3,4), base_modes(3,5));
fprintf('  Mode 3: f = %.4f MHz, Q = %.0f, R/Q = %.0f Ohm/m^2\n', ...
    base_modes(4,3)/1e6, base_modes(4,4), base_modes(4,5));

%% ========== RUN BBU TRACKING ==========
p = struct();
p.machine_index     = 0;
p.active_modes      = [2, 3];      % 0-based mode indices
p.f0                = 1.3e9;       % [Hz]
p.Nc                = Nc;
p.Na                = Na;
p.NmTotal           = NmTotal;
p.Nb                = 20000;       % number of bunches (increase for better statistics)
p.I_beam            = 1.0;         % beam current [A] — try different values
p.initial_V_dip     = 1.0;         % HOM seed voltage [V]
p.long_motion       = false;
p.output_dir        = 'results';
p.max_volt_lines    = 200;
p.output_time_units = 1e-6;        % output in microseconds
p.machines_dir      = machines_dir;

fprintf('\n--- Running BBU dipole tracking ---\n');
fprintf('  Nb = %d,  I = %.3f A\n\n', p.Nb, p.I_beam);

results = BBU_dipole_tracking(p);

%% ========== POST-PROCESSING ==========
fprintf('\n--- Summary ---\n');
if results.stable
    fprintf('  STABLE at I = %.3f A\n', p.I_beam);
    fprintf('  All HOM voltages are damped.\n');
else
    fprintf('  UNSTABLE at I = %.3f A\n', p.I_beam);
    fprintf('  At least one HOM voltage is growing.\n');
end

fprintf('\n  Growth rates:\n');
for ic = 1:Nc
    for im = 1:length(p.active_modes)
        fprintf('    Cav %d, Mode %d: decrement = %+.4e 1/s\n', ...
            ic, p.active_modes(im), results.decrement(ic, im));
    end
end

fprintf('\nTo scan the threshold current, try different p.I_beam values.\n');


%% ========== HELPER: FODO transport matrix ==========
function M = fodo_matrix(mu, L)
% FODO_MATRIX  6x6 transport matrix for a FODO cell
%   mu : betatron phase advance [rad]
%   L  : total cell length [m]
% Returns a block-diagonal matrix with equal x and y focusing.

    cs = cos(mu);
    sn = sin(mu);
    M = eye(6);
    M(1,1) = cs;               M(1,2) = L * sn / mu;
    M(2,1) = -mu * sn / L;     M(2,2) = cs;
    M(3,3) = cs;               M(3,4) = L * sn / mu;
    M(4,3) = -mu * sn / L;     M(4,4) = cs;
end
