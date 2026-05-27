function Z = csr_drift(kC, s_p, dip, egamma)
% csr_drift  CSR drift impedance 
%
% After the beam exits a dipole, the coherent radiation wake continues
% to act on the beam in the downstream drift.
% Ref: R. Bosch, PRST-AB 11, 090702 (2008)
%
% Inputs:
%   kC     - compressed wavenumber k*C(s') [m^-1]
%   s_p    - position [m] (must be in drift, outside dipole)
%   dip    - dipole table struct (.N, .s1, .s2, .rho) all [m]
%   egamma - Lorentz gamma

    if dip.N == 0, Z = 0; return; end

    % Collect all upstream dipoles
    Lb  = [];  rho = [];  ds = [];
    for m = 1:dip.N
        if s_p > dip.s2(m)
            Lb(end+1)  = abs(dip.s2(m) - dip.s1(m));   %#ok
            rho(end+1) = dip.rho(m);                    %#ok
            ds(end+1)  = s_p - dip.s2(m);               %#ok
        end
    end
    if isempty(Lb), Z = 0; return; end

    Lb = Lb(:);  rho = rho(:);  ds = ds(:);
    r  = length(ds);

    % Case D: Bosch Eq.(24), PRST-AB 10, 050701 (2007)
    tmp03 = zeros(r, 1);
    for m = 1:r
        if ds(m) > abs(rho(m))^(2/3) * (2*pi/kC)^(1/3)
            d_eff    = min(ds(m), egamma^2/kC);
            tmp03(m) = 2 / d_eff;
        end
    end

    % Case C
    tmp04 = exp(-1i*kC*Lb.^2 ./ (6*rho.^2) .* (Lb + 3*ds));
    tmp06 = -4*tmp04 ./ (Lb + 2*ds);

    Z = tmp03(end) + tmp06(end);
end
