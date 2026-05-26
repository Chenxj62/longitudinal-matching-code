function GenerateMatrices(beamline, E0_eV, varargin)
% GenerateMatrices  Compute BBU inter-cavity transport matrices from a beamline
%
% Replaces GenerateMachines.py (Elegant-based) in the UniBBU workflow.
% Finds all cavity passages in a beamline, computes the 6x6 transport
% matrix and arc length between consecutive cavities, and writes
% machines/matrices/0000.txt for BBU_dipole_tracking.m.
%
% Quad binding (same as TSCenv.m):
%   G_norm = k1 * gamma  is energy-invariant for a fixed physical magnet.
%   Quads with the SAME name share the same G_norm (set at first pass).
%   At each arc, k1_eff = G_norm / gamma_current.
%
% Units (consistent with createElement.m / getElementMatrix.m):
%   gradient : MV/m  ->  energy gain = gradient * L * cos(phi)  [MeV]
%   phase    : degrees
%   length   : m,  k1 : 1/m^2
%
% -----------------------------------------------------------------------
% USAGE
%   GenerateMatrices(beamline, E0_eV)
%   GenerateMatrices(beamline, E0_eV, 'MachineIndex', 0)
%   GenerateMatrices(beamline, E0_eV, 'OutputDir', 'machines/matrices')
%
% Required:
%   beamline  - struct array (from createElement / readBmad).
%               Must be the FULL sequence for all passes in order.
%   E0_eV     - beam energy at beamline entrance [eV]
%
% Optional name-value:
%   'MachineIndex'  integer   file index     (default 0  ->  0000.txt)
%   'OutputDir'     string    output folder  (default 'machines/matrices')
%   'CavityTypes'   cell      type strings   (default {'cavity','lcavity','rf','rfcavity'})
%   'Verbose'       logical   print table    (default true)
%
% -----------------------------------------------------------------------
% MULTI-PASS ERL EXAMPLE (Nc=2 cavities, 4 linac passes, Na=8 passages)
%
%   linac = [d, q, cav1, d, cav2, q, d];   % built with createElement
%   arc1  = [d, q, d, q, d];               % recirculation arc 1
%   arc2  = [d, q, d, q, d];               % recirculation arc 2
%   arcR  = [d, q, d];                     % return arc
%
%   full = [linac, arc1, linac, arc2, linac, arcR, linac];
%   GenerateMatrices(full, 5e6);
%
% -----------------------------------------------------------------------
% OUTPUT FILE FORMAT  (machines/matrices/0000.txt)
%   For each arc ia = 1..Na:
%     Line 1   : Ea [eV]  beam energy at entrance of arc ia
%     Line 2   : La [m]   arc length
%     Lines 3-8: 6x6 transport matrix, one row per line, tab-separated
% -----------------------------------------------------------------------

%% ===== OPTIONS =====
opt = inputParser;
addRequired (opt, 'beamline');
addRequired (opt, 'E0_eV',       @isnumeric);
addParameter(opt, 'MachineIndex', 0,                   @isnumeric);
addParameter(opt, 'OutputDir',   'machines/matrices',  @ischar);
addParameter(opt, 'CavityTypes', {'cavity','lcavity','rf','rfcavity'}, @iscell);
addParameter(opt, 'Verbose',     true,                 @islogical);
parse(opt, beamline, E0_eV, varargin{:});

machine_idx = opt.Results.MachineIndex;
output_dir  = opt.Results.OutputDir;
cav_types   = lower(opt.Results.CavityTypes);
verbose     = opt.Results.Verbose;

me_eV = 0.511e6;   % electron rest energy [eV]

%% ===== PRE-PASS: QUAD BINDING =====
% Same logic as TSCenv.m:
%   G_norm = k1 * gamma  (energy-invariant for a fixed physical magnet)
% First occurrence of each quad name registers G_norm.
% Repeated occurrences at different energies get k1 = G_norm / gamma.

quad_name_map = containers.Map('KeyType','char','ValueType','double');
temp_E = E0_eV;

for i = 1:length(beamline)
    elem  = beamline(i);
    gamma = temp_E / me_eV;

    if ismember(lower(elem.type), cav_types)
        grad = 0; if isfield(elem,'gradient'), grad = elem.gradient; end
        phi  = 0; if isfield(elem,'phase'),    phi  = elem.phase * pi/180; end
        temp_E = temp_E + grad * elem.length * cos(phi) * 1e6;   % eV

    elseif strcmpi(elem.type, 'quadrupole')
        qname = '';
        if isfield(elem,'name'), qname = elem.name; end
        if ~isempty(qname)
            G_norm = elem.k1 * gamma;
            if quad_name_map.isKey(qname)
                beamline(i).k1 = quad_name_map(qname) / gamma;
            else
                quad_name_map(qname) = G_norm;
            end
        end
    end
end

if verbose
    fprintf('Quad binding: %d unique quad names.\n', quad_name_map.Count);
end

%% ===== SCAN: SPLIT BEAMLINE INTO ARCS =====
% arc_elements{ia} = beamline elements BEFORE cavity ia
% cav_L, cav_grad, cav_phi : cavity parameters collected during scan

arc_elements = {};
cav_L        = [];
cav_grad_arr = [];
cav_phi_arr  = [];
seg_start    = 1;

for k = 1:length(beamline)
    if ismember(lower(beamline(k).type), cav_types)

        % Close current arc
        if k > seg_start
            arc_elements{end+1} = beamline(seg_start : k-1); %#ok<AGROW>
        else
            arc_elements{end+1} = createElement('drift', 0);  %#ok<AGROW>
        end

        % Record cavity kick parameters
        g   = 0; if isfield(beamline(k),'gradient'), g   = beamline(k).gradient;          end
        phi = 0; if isfield(beamline(k),'phase'),    phi = beamline(k).phase * pi / 180;   end
        cav_L(end+1)        = beamline(k).length;   %#ok<AGROW>
        cav_grad_arr(end+1) = g;                    %#ok<AGROW>
        cav_phi_arr(end+1)  = phi;                  %#ok<AGROW>

        seg_start = k + 1;
    end
end

Na = length(arc_elements);
if Na == 0
    types = strjoin(unique({beamline.type}), ', ');
    error('No cavities found. Types present: %s\nAdjust ''CavityTypes'' if needed.', types);
end

%% ===== MAIN PASS: MATRICES AND ENERGY TRACKING =====
Ea = zeros(Na, 1);
La = zeros(Na, 1);
Ta = zeros(6, 6, Na);

E = E0_eV;

for ia = 1:Na
    gamma = E / me_eV;
    seg   = arc_elements{ia};

    % Apply quad binding to this arc's elements
    for j = 1:length(seg)
        if strcmpi(seg(j).type, 'quadrupole')
            qname = '';
            if isfield(seg(j),'name'), qname = seg(j).name; end
            if ~isempty(qname) && quad_name_map.isKey(qname)
                seg(j).k1 = quad_name_map(qname) / gamma;
            end
        end
    end

    % Transport matrix for this arc (gamma constant, no cavities inside)
    [R_arc, ~] = calculateTotalMatrix(seg, gamma);

    Ea(ia)     = E;
    La(ia)     = sum([seg.length]);
    Ta(:,:,ia) = R_arc;

    % Energy kick from cavity ia
    dE_eV = cav_grad_arr(ia) * cav_L(ia) * cos(cav_phi_arr(ia)) * 1e6;
    E     = E + dE_eV;

    if E < me_eV
        warning('Arc %d: energy %.3e eV < rest mass after cavity. Check phase/gradient.', ia, E);
        E = me_eV;
    end
end

%% ===== PRINT SUMMARY =====
if verbose
    fprintf('\n  Arc   E_in [MeV]   L [m]   Grad [MV/m]   phi [deg]   dE [MeV]   E_out [MeV]\n');
    fprintf(  '  ---   ----------   -----   -----------   ---------   --------   -----------\n');
    E_run = E0_eV;
    for ia = 1:Na
        dE    = cav_grad_arr(ia) * cav_L(ia) * cos(cav_phi_arr(ia));
        fprintf('  %3d   %10.4f   %5.3f   %11.4f   %9.2f   %8.4f   %11.4f\n', ...
            ia, E_run/1e6, La(ia), cav_grad_arr(ia), cav_phi_arr(ia)*180/pi, dE, E_run/1e6+dE);
        E_run = E_run + dE * 1e6;
    end
    fprintf('\n');
end

%% ===== WRITE OUTPUT FILE =====
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
outfile = fullfile(output_dir, sprintf('%04d.txt', machine_idx));

fid = fopen(outfile, 'w');
if fid < 0, error('Cannot open for writing: %s', outfile); end

for ia = 1:Na
    fprintf(fid, '%e\n', Ea(ia));
    fprintf(fid, '%e\n', La(ia));
    for row = 1:6
        fprintf(fid, '%e\t', Ta(row,:,ia));
        fprintf(fid, '\n');
    end
end
fclose(fid);

if verbose, fprintf('Written: %s   (%d arcs)\n', outfile, Na); end

end
