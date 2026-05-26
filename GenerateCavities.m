function GenerateCavities(hom_file, Nc, varargin)
% GenerateCavities  Convert a nominal HOM table into machines/cavities/XXXX.txt
%
% Reads a HOM parameter file (e.g. TESLA1.3GHz.txt) and writes the cavity
% parameter file required by BBU_dipole_tracking.m.
%
% For the time-domain tracking, all Nc cavities receive the same nominal
% HOM parameters (no randomization). Optional Gaussian frequency spread
% and Gamma Q-factor spread can be enabled if statistical studies are needed.
%
% -----------------------------------------------------------------------
% INPUT FILE FORMAT  (e.g. TESLA1.3GHz.txt)  — 8 columns, space/tab separated:
%   index   order   freq[Hz]   Q   R/Q[Ohm/m^2]   theta[deg]   df[Hz]   dQQ
%
% OUTPUT FILE FORMAT  (machines/cavities/0000.txt)  — 6 columns per line:
%   index   order   freq[Hz]   Q   R/Q[Ohm/m^2]   theta[deg]
%   (Nc cavities × NmTotal modes, written consecutively)
%
% -----------------------------------------------------------------------
% USAGE
%   GenerateCavities('TESLA1.3GHz.txt', 4)
%   GenerateCavities('TESLA1.3GHz.txt', 4, 'MachineIndex', 0)
%   GenerateCavities('TESLA1.3GHz.txt', 4, 'FreqSpread', 1e6, 'NumMachines', 10)
%
% Required:
%   hom_file  - path to nominal HOM parameter file
%   Nc        - number of cavities
%
% Optional name-value:
%   'MachineIndex'   integer   single machine index to write  (default 0)
%   'NumMachines'    integer   write multiple randomized machines (default 1)
%   'OutputDir'      string    output folder  (default 'machines/cavities')
%   'FreqSpread'     double    1-sigma frequency spread [Hz] (default 0 = nominal)
%   'QSpread'        double    relative Q spread, sigma of Gamma dist (default 0 = nominal)
%   'CutoffSigma'    double    truncation of spread distributions (default 3)
%   'Verbose'        logical   print summary (default true)
% -----------------------------------------------------------------------

%% ===== OPTIONS =====
opt = inputParser;
addRequired (opt, 'hom_file', @ischar);
addRequired (opt, 'Nc',       @isnumeric);
addParameter(opt, 'MachineIndex', 0,                  @isnumeric);
addParameter(opt, 'NumMachines',  1,                  @isnumeric);
addParameter(opt, 'OutputDir',   'machines/cavities', @ischar);
addParameter(opt, 'FreqSpread',   0,                  @isnumeric);
addParameter(opt, 'QSpread',      0,                  @isnumeric);
addParameter(opt, 'CutoffSigma',  3,                  @isnumeric);
addParameter(opt, 'Verbose',      true,               @islogical);
parse(opt, hom_file, Nc, varargin{:});

start_idx   = opt.Results.MachineIndex;
num_machines = opt.Results.NumMachines;
output_dir  = opt.Results.OutputDir;
df_sigma    = opt.Results.FreqSpread;
dQQ_sigma   = opt.Results.QSpread;
cutoff      = opt.Results.CutoffSigma;
verbose     = opt.Results.Verbose;

%% ===== READ NOMINAL HOM FILE =====
fid = fopen(hom_file, 'r');
if fid < 0, error('Cannot open HOM file: %s', hom_file); end

raw = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || line(1) == '!' || line(1) == '%'
        continue;
    end
    vals = sscanf(line, '%f');
    if length(vals) >= 6
        raw(end+1, :) = vals(1:8)'; %#ok<AGROW>  % keep all 8 cols
    end
end
fclose(fid);

if isempty(raw)
    error('No data read from %s. Check file format.', hom_file);
end

NmTotal = size(raw, 1);   % total modes per cavity

% Nominal values
n_idx  = raw(:,1);   % mode index
order  = raw(:,2);   % azimuthal order
f0     = raw(:,3);   % frequency [Hz]
Q0     = raw(:,4);   % Q factor
RoQ    = raw(:,5);   % R/Q [Ohm/m^2]
theta0 = raw(:,6);   % polarization angle [deg]
df_nom = raw(:,7);   % frequency spread (from file, used if FreqSpread not set)
dQQ_nom = raw(:,8);  % Q spread (from file)

% Use file spread values if not overridden
if df_sigma == 0 && any(df_nom ~= 0)
    df_sigma = 0;    % keep nominal (no randomization) unless user asks
end

if verbose
    fprintf('Read %d modes from %s\n', NmTotal, hom_file);
    orders = unique(order)';
    for o = orders
        fprintf('  Order %d: %d modes\n', o, sum(order == o));
    end
end

%% ===== WRITE CAVITY FILES =====
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

for mi = start_idx : start_idx + num_machines - 1

    outfile = fullfile(output_dir, sprintf('%04d.txt', mi));
    fid = fopen(outfile, 'w');
    if fid < 0, error('Cannot write: %s', outfile); end

    for ic = 1:Nc
        for k = 1:NmTotal

            % --- Frequency: nominal + optional Gaussian spread ---
            if df_sigma > 0
                df = sampleTruncatedGaussian(0, df_sigma, cutoff);
            else
                df = 0;
            end
            f_k = f0(k) + df;

            % --- Q factor: nominal × optional Gamma-distributed factor ---
            if dQQ_sigma > 0
                % Gamma(2, 0.5) has mean=1, truncated at 1 ± cutoff/sqrt(2)
                dQQ = sampleTruncatedGamma(cutoff);
            else
                dQQ = 1;
            end
            Q_k = Q0(k) * dQQ;

            % --- Polarization: keep nominal ---
            theta_k = theta0(k);

            fprintf(fid, '%d\t%d\t%e\t%e\t%e\t%e\n', ...
                n_idx(k), order(k), f_k, Q_k, RoQ(k), theta_k);
        end
    end

    fclose(fid);
    if verbose, fprintf('Written: %s  (%d cavities × %d modes)\n', outfile, Nc, NmTotal); end
end

end  % GenerateCavities

%% =========================================================================
function x = sampleTruncatedGaussian(mu, sigma, cutoff)
% Sample from a Gaussian truncated at ±cutoff*sigma
x = cutoff * sigma + 1;
while abs(x - mu) > cutoff * sigma
    x = mu + sigma * randn();
end
end

%% =========================================================================
function dQQ = sampleTruncatedGamma(cutoff)
% Sample Gamma(2, 0.5) truncated at 1 ± cutoff/sqrt(2)  (same as GenerateMachines.py)
limit = 1 + cutoff / sqrt(2);
dQQ = limit + 1;
while abs(dQQ) > limit
    dQQ = gamrnd(2, 0.5);
end
end
