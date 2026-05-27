function [z_out, delta_out, E_ref_out] = cavity_long(z_in, delta_in, ...
                                                     E_ref_in, V, phi_deg, f, ...
                                                     series_order)
% CAVITY_LONG  Thin-cavity longitudinal evolution (slice version).
%
%   [z_out, delta_out, E_ref_out] = ...
%       cavity_long(z_in, delta_in, E_ref_in, V, phi_deg, f, series_order)
%
%   Applies an RF cavity to a slice-binned beam.  Inputs are slice arrays
%   (slice center z, slice-mean delta); the formula is element-wise so
%   the slice current I_slice is unchanged and need not be passed in or out.
%
%   Inputs
%     z_in         [Ns x 1]   slice center z [m]
%     delta_in     [Ns x 1]   slice-mean relative energy deviation
%     E_ref_in     scalar     reference energy entering the cavity [MeV]
%     V            scalar     cavity voltage [MeV]
%     phi_deg      scalar     RF phase [deg]  (on-crest = 0)
%     f            scalar     RF frequency [Hz]
%     series_order scalar     (optional) Taylor-expand cos(phi-kz)-cos(phi)
%                             in u = k*z to this order. 0 = use full cosine
%                             (default). Typical values: 2 (chirp+curvature),
%                             3, 4.
%
%   Outputs
%     z_out      = z_in                                   (thin cavity)
%     delta_out  cosine RF kick (full or Taylor-truncated) + adiabatic damping
%     E_ref_out  E_ref_in + V*cos(phi)                    [MeV]
%
%   Formula (series_order = 0):
%       E_out     = E_in + V*cos(phi)
%       delta_out = (E_in / E_out) * delta_in
%                 + (V    / E_out) * (cos(phi - k*z_in) - cos(phi))
%       z_out     = z_in
%   with k = 2*pi*f/c. Phase convention: phi = 0 -> on-crest.
%
%   Series expansion (series_order = N >= 1):
%       cos(phi - u) - cos(phi)
%         = sum_{n=1..N} c_n / n! * u^n,  u = k*z
%       c_n cycle (n mod 4): 1->sin(phi), 2->-cos(phi), 3->-sin(phi), 0->cos(phi)
%
%   See also: ../mp/cavity_long.m (per-particle version).

    if nargin < 7 || isempty(series_order)
        series_order = 0;
    end

    c       = 2.99792458e8;
    k       = 2*pi * f / c;
    phi_rad = deg2rad(phi_deg);

    E_ref_out = E_ref_in + V * cos(phi_rad);

    if series_order == 0
        kick = cos(phi_rad - k*z_in) - cos(phi_rad);
    else
        u    = k * z_in;
        s    = sin(phi_rad);
        co   = cos(phi_rad);
        kick = zeros(size(u));
        un   = ones(size(u));        % u^0
        for n = 1:series_order
            un = un .* u;            % u^n
            switch mod(n, 4)
                case 1, c_n =  s;
                case 2, c_n = -co;
                case 3, c_n = -s;
                case 0, c_n =  co;
            end
            kick = kick + (c_n / factorial(n)) .* un;
        end
    end

    delta_out = (E_ref_in / E_ref_out) .* delta_in ...
              + (V        / E_ref_out) .* kick;

    z_out = z_in;
end
