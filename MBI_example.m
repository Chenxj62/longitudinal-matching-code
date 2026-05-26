%% MBI_example - How to use MBI_calc with your beamline
% Uses your createElement / getElementMatrix format directly.
% No ELEGANT needed.

%% Step 1: Build your beamline with createElement
% Example: a simple bunch compressor (4-dipole chicane)

% Drifts
d1 = createElement('drift', 1.0);
d2 = createElement('drift', 0.5);
d3 = createElement('drift', 2.0);

% Dipoles (chicane: +angle, -angle, -angle, +angle)
rho = 5.0;           % bending radius [m]
L_bend = 0.5;        % dipole length [m]
h = 1/rho;            % curvature [m^-1]

b1 = createElement('dipole', L_bend, h);
b2 = createElement('dipole', L_bend, -h);
b3 = createElement('dipole', L_bend, -h);
b4 = createElement('dipole', L_bend, h);

% Quadrupoles
q1 = createElement('quadrupole', 0.2, 1.5);    % k1 = 1.5 m^-2
q2 = createElement('quadrupole', 0.2, -1.5);

% Assemble beamline
beamline = [d1, q1, d2, b1, d3, b2, d2, d2, b3, d3, b4, d2, q2, d1];

%% Step 2: Set beam parameters
beam.energy_MeV  = 250;       % beam energy [MeV]
beam.I_b         = 30;        % peak current [A]
beam.emit_norm_x = 1.0;       % normalized emittance [um]
beam.emit_norm_y = 1.0;       % normalized emittance [um]
beam.betax0      = 10.0;      % initial betax [m]
beam.betay0      = 10.0;      % initial betay [m]
beam.alphax0     = 0.0;
beam.alphay0     = 0.0;
beam.sigma_delta = 1e-4;      % relative energy spread
beam.chirp       = 0.0;       % chirp [m^-1]
beam.sigma_z     = 0.5e-3;    % rms bunch length [m]

%% Step 3: Set options
opts.lambda_range     = [1, 200];   % wavelength range [um]
opts.scan_num         = 50;
opts.mesh_num         = 400;
opts.iCSR_ss          = 1;          % steady-state CSR on
opts.iCSR_tr          = 0;
opts.iCSR_drift       = 0;
opts.iLSC             = 0;
opts.calc_energy_gain = 1;
opts.plot             = 1;

%% Step 4: Run
result = MBI_calc(beamline, beam, opts);

%% Step 5: Use results
fprintf('\n=== Results ===\n');
fprintf('Peak density gain: %.4f at lambda = %.2f um\n', ...
    max(result.Gf_dd), result.lambda(result.Gf_dd == max(result.Gf_dd)));

% Available fields:
%   result.lambda     - wavelength [um]
%   result.s          - position [m]
%   result.Gf_dd      - density-to-density gain spectrum
%   result.Gf_de      - density-to-energy gain spectrum
%   result.Gf_ed      - energy-to-density gain spectrum
%   result.Gf_total   - total density gain
%   result.G_vs_s     - gain map G(s, lambda)
