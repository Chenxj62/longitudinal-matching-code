function Z = lsc_impedance(k, sig_x, sig_y, egamma, model)
% lsc_impedance  1D longitudinal space charge impedance (SI units)
%
% Inputs:
%   k      - wavenumber [m^-1] (already includes compression: k*C)
%   sig_x  - rms horizontal beam size [m]
%   sig_y  - rms vertical beam size [m]
%   egamma - Lorentz gamma
%   model  - 1: on-axis uniform circular
%            2: averaged uniform circular
%            3: transverse Gaussian (default)
%
% Output:
%   Z      - LSC impedance [m^-1] (same units as CSR impedances)
%
% Reference: Huang, Stupakov, PRST-AB 7, 015001 (2004)
%            Venturini, PRST-AB 11, 034401 (2008)

    if nargin < 5, model = 3; end

    rb = 1.747/2 * (sig_x + sig_y);   % effective beam radius [m]
    xi = k * rb / egamma;

    switch model
        case 1
            % On-axis, transverse uniform circular cross section
            Z = 4i / (rb * egamma) * (1 - xi*besselk(1,xi)) / xi;

        case 2
            % Averaged over transverse distribution, uniform circular
            Z = 4i / (rb * egamma) * (1 - 2*besseli(1,xi)*besselk(1,xi)) / xi;

        case 3
            % Transverse axisymmetric Gaussian
            S_avg  = 0.5 * (sig_x + sig_y);
            xi_a   = k * S_avg / egamma;
            tmp01  = xi_a / (egamma * S_avg);
            tmp02  = 0.5 * xi_a^2;
            tmp03  = exp(tmp02);
            tmp04  = -expint(tmp02);   % Matlab expint = E1(z)
            Z = -1i * tmp01 * tmp03 * tmp04;

        otherwise
            error('lsc_impedance: unknown model %d', model);
    end
end
