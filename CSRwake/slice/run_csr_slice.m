% RUN_CSR_SLICE  Example driver for csr_calc_slice.m
%
%   Workflow:
%     1. readBmad(latticefile)   ->  beamline (struct array of elements)
%     2. bmad2slice(bunchfile)   ->  (z_slice, delta_slice, I_slice)
%     3. csr_calc_slice(...)     ->  end-of-line phase space + currents
%
%   `getElementMatrix.m`, `readBmad.m`, `bmad2slice.m`, and `createElement.m`
%   live in the parent code directory; we addpath them at the top.

clear; clc; close all;

addpath('..');         % bmad2slice.m
addpath('..\..');      % readBmad.m, getElementMatrix.m, createElement.m

% =================== USER INPUTS ===================
latticefile = '..\..\latv13.arc';        % .bmad lattice file
bunchfile   = '..\Bunch.ast';           % ASTRA particle file
Q_total     = 70e-12;                    % [C] total bunch charge

% slice extraction (bmad2slice)
nslices       = 200;
fit_order     = 7;
min_particles = 10;

% CSR / transport
gamma_rel = 1000 / 0.511;                % gamma at this arc's entrance
T566      = -1.1;                        % [m] 2nd-order longitudinal dispersion

% smoothing inside csr_calc_slice
smooth_win = 50;

% override sigzf if you want to bypass the auto value (empty = auto)
sigzf_override = 12e-6;

% optional dR56 added only to the transport step (does NOT enter S integrals)
dR56 = 0;



% R56(s) substeps per dipole inside the line integral
n_substeps = 200;
% ===================================================

% --- 1) lattice from BMAD file ----------------------------------------
beamline = readBmad(latticefile);

% --- 2) slice extraction from the bunch ASTRA file --------------------
%   bmad2slice returns: z_slice, delta_slice (fitted), I_slice
%   Use source_file = Q_total (in C) for the ASTRA-style branch.
[z_slice, delta_slice, I_slice] = bmad2slice( ...
    bunchfile, nslices, fit_order, min_particles, Q_total);

% --- 3) CSR + transport on the slice arrays ---------------------------
out = csr_calc_slice(z_slice, delta_slice, I_slice, ...
        beamline, gamma_rel, T566, ...
        'smooth_win', smooth_win, ...
        'dR56',       dR56, ...
        'n_substeps', n_substeps, ...
        'sigzf',      sigzf_override, ...
        'plot_flag',  1);

% --- 4) print resolved scalars ----------------------------------------
ud = out.Properties.UserData;
fprintf(['Resolved scalars:\n' ...
         '  sigma_z0 = %.3g um   sigma_zf = %.3g um\n' ...
         '  sige     = %.3g     r56_tot   = %.3g m\n' ...
         '  S_R56    = %.3g     S_one     = %.3g\n' ...
         '  n_bends_used = %d\n'], ...
        ud.sigma_z0*1e6, ud.sigma_zf*1e6, ud.sige_used, ud.r56_tot, ...
        ud.S_R56, ud.S_one, ud.n_bends_used);
