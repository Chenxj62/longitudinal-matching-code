function Z = csr_ss(kC, rho, A_csr)
% csr_ss  Steady-state CSR impedance (ultra-relativistic, free-space)
%
%   Z = -i * (kC)^{1/3} * A * |rho|^{-2/3}
%
% Inputs:
%   kC    - compressed wavenumber k*C(s') [m^-1]
%   rho   - bending radius [m]
%   A_csr - CSR constant = 3^{-1/3} * Gamma(2/3) * (i*sqrt(3) - 1)

    Z = -1i * kC^(1/3) * A_csr * abs(rho)^(-2/3);
end
