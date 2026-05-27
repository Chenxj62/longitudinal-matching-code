function Z = compute_Z_csr(k, C_sp, s_p, rho_p, dip, egamma, opts, A_csr)
% compute_Z_csr  Total CSR impedance Z(k*C(s'), s') — SI units [m]
%
% Inputs:
%   k      - wavenumber [m^-1]
%   C_sp   - compression factor at s'
%   s_p    - position s' [m]
%   rho_p  - bending radius at s' [m] (1e50 if drift)
%   dip    - dipole table struct (.N, .s1, .s2, .rho) all [m]
%   egamma - Lorentz gamma
%   opts   - option struct
%   A_csr  - CSR constant

    kC = k * C_sp;
    Z  = 0;
    in_dipole = abs(rho_p) < 1e20;

    % Steady-state CSR
    if opts.iCSR_ss && in_dipole
        if opts.issCSRpp
            h_pipe = opts.full_pipe_height;        % [m]
            k_th   = 2*pi * sqrt(abs(rho_p) / h_pipe^3);
            if kC < k_th
                Z = Z + csr_shielding(kC, h_pipe, rho_p);
            else
                Z = Z + csr_ss(kC, rho_p, A_csr);
            end
        else
            Z = Z + csr_ss(kC, rho_p, A_csr);
        end
    end

    % Transient CSR
    if opts.iCSR_tr && in_dipole
        Z = Z + csr_tr(kC, s_p, rho_p, dip);
    end

    % CSR drift
    if opts.iCSR_drift && ~in_dipole
        Z = Z + csr_drift(kC, s_p, dip, egamma);
    end
end
