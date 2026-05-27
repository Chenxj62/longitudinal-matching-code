%% MBI_example - How to use MBI_calc (SI units)
%  All inputs in SI: meters, radians, amperes.
close all;
clc;

%% Add parent path (where createElement, getElementMatrix live)
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

beamline = readBmad("lata3v12.bmad");

%% Step 2: Beam parameters (all SI)
beam.energy_MeV  = 1000;
beam.I_b         = 7;            % [A]
beam.emit_norm_x = 1.01e-6;      % [m-rad]  (1.01 um)
beam.emit_norm_y = 1.01e-6;      % [m-rad]
beam.betax0      = 45.8;         % [m]
beam.betay0      = 140.9;        % [m]
beam.alphax0     = 2.37;
beam.alphay0     = 2.89;
beam.sigma_delta = 1.157e-6;
beam.chirp       = 3.33;          % [m^-1]
beam.sigma_z     = 0.648e-3;     % [m]
initial_beta=[beam.betax0,beam.betay0,beam.alphax0,beam.alphay0];
%calcBetaEvolution(beamline, 1000, initial_beta, 1)
%% Step 3: Options (lambda_range in meters)
opts.lambda_range     = [2e-6, 2e-3];  % [m]  (2 um - 2 mm)
opts.scan_num         = 500;
opts.mesh_num         = 400;
opts.iCSR_ss          = 1;
opts.iCSR_tr          = 1;
opts.iCSR_drift       = 1;
opts.iLSC             = 0;
opts.LSC_model        = 0;
opts.calc_energy_gain = 1;
opts.plot             = 1;

%% Step 4: Run
result = MBI_calc(beamline, beam, opts);

%% Step 5: Results
fprintf('\n=== Results ===\n');
[pk, idx] = max(result.Gf_dd);
fprintf('Peak density gain: %.4f at lambda = %.4e m (%.2f um)\n', ...
    pk, result.lambda(idx), result.lambda(idx)*1e6);
