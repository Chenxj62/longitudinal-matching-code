function [z_out, delta_out, E_ref_out] = cavity_long(z_in, delta_in, ...
                                                     E_ref_in, V, phi_deg, f)
% CAVITY_LONG  Thin-cavity longitudinal evolution (multi-particle version).
%
%   [z_out, delta_out, E_ref_out] = ...
%       cavity_long(z_in, delta_in, E_ref_in, V, phi_deg, f)
%
%   Applies an RF cavity to a full particle distribution.  Inputs are
%   per-particle arrays of length N; the formula is element-wise.
%
%   Inputs
%     z_in       [N x 1]    particle longitudinal coordinate [m]
%     delta_in   [N x 1]    particle relative energy deviation
%     E_ref_in   scalar     reference energy entering the cavity [MeV]
%     V          scalar     cavity voltage [MeV]
%     phi_deg    scalar     RF phase [deg]  (on-crest = 0)
%     f          scalar     RF frequency [Hz]
%
%   Outputs
%     z_out      = z_in                                   (thin cavity)
%     delta_out  full nonlinear cosine RF kick with adiabatic damping
%     E_ref_out  E_ref_in + V*cos(phi)                    [MeV]
%
%   Formula:
%       E_out     = E_in + V*cos(phi)
%       delta_out = (E_in / E_out) * delta_in
%                 + (V    / E_out) * (cos(phi - k*z_in) - cos(phi))
%       z_out     = z_in
%   with k = 2*pi*f/c. Phase convention: phi = 0 -> on-crest.
%
%   See also: ../slice/cavity_long.m (slice version).

    c       = 2.99792458e8;
    k       = 2*pi * f / c;
    phi_rad = deg2rad(phi_deg);

    E_ref_out = E_ref_in + V * cos(phi_rad);

    delta_out = (E_ref_in / E_ref_out) .* delta_in ...
              + (V        / E_ref_out) .* (cos(phi_rad - k*z_in) - cos(phi_rad));

    z_out = z_in;
end
