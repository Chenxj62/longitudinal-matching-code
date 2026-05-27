function Z = csr_shielding(kC, h_pipe, rho)
% csr_shielding  CSR impedance with parallel-plate shielding
%
% At long wavelengths (k < k_th), the conducting pipe walls screen
% the CSR field. Uses the Agoh-Yokoya parallel plate model.
% Ref: Agoh & Yokoya, PRSTAB 7, 054403 (2004)
%
% Inputs:
%   kC     - compressed wavenumber k*C(s') [m^-1]
%   h_pipe - full pipe height [m]
%   rho    - bending radius [m]

    p_max   = 200;
    rho_abs = abs(rho);
    coeff   = 8*pi^2/h_pipe * (2/(kC*rho_abs))^(1/3);
    s = 0;
    for p = 0:p_max
        beta_p = (2*p+1)*pi/h_pipe * (rho_abs/(2*kC^2))^(1/3);
        s = s + F0_shielding(beta_p);
    end
    Z = coeff * s;
end
