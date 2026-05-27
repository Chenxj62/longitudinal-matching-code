function Z = csr_tr(kC, s_p, rho, dip)
% csr_tr  Transient CSR impedance (entrance transient inside dipole)
%
% Within formation length L_f = (24*lambda*rho^2)^{1/3}, the steady-state
% CSR has not yet developed. This function returns the transient impedance.
%
% Inputs:
%   kC   - compressed wavenumber k*C(s') [m^-1]
%   s_p  - position [m]
%   rho  - bending radius [m]
%   dip  - dipole table struct (.N, .s1, .s2, .rho) all [m]

    local_s = local_s_in_dipole(s_p, dip);
    if local_s == 0,    Z = 0; return; end
    if local_s < 5e-4,  local_s = 5e-4; end    % 0.5 mm minimum

    lam_tmp = 2*pi / kC;
    L_form  = (24 * lam_tmp * rho^2)^(1/3);

    if local_s < L_form
        z_L = local_s^3 / (24 * rho^2);
        mu  = kC * z_L;
        t1  = -4/local_s * exp(-1i*4*mu);
        t2  = 4/(3*local_s) * (1i*mu)^(1/3) * upper_igamma(-1/3, 1i*mu);
        Z   = t1 + t2;
    else
        Z = 0;
    end
end


% Upper incomplete gamma: Gamma(a,z) for complex z
% Replaces mfun('GAMMA', a, z) without Symbolic Math Toolbox
function val = upper_igamma(a, z)
    ga = exp(gammalnz_local(a));
    val = ga * (1 - gammaincc_local(z, a));
end

function b = gammaincc_local(x, a)
    if x == 0, b = 0; return; end
    lnga = gammalnz_local(a);
    if abs(x) < abs(a)+1
        ap = a; s = 1/a; del = s;
        for n = 1:500
            ap = ap + 1;
            del = x * del / ap;
            s = s + del;
            if abs(del) < 10*eps*abs(s), break; end
        end
        b = s * exp(-x + a*log(x) - lnga);
    else
        a0 = 1; a1 = x; b0 = 0; b1 = 1;
        fac = 1; g = b1; gold = b0; n = 1;
        while abs(g - gold) >= 10*eps*abs(g)
            gold = g;
            ana = n - a;
            a0 = (a1 + a0*ana) * fac;
            b0 = (b1 + b0*ana) * fac;
            anf = n * fac;
            a1 = x * a0 + anf * a1;
            b1 = x * b0 + anf * b1;
            fac = 1 / a1;
            g = b1 * fac;
            n = n + 1;
            if n > 500, break; end
        end
        b = 1 - exp(-x + a*log(x) - lnga) * g;
    end
end

function f = gammalnz_local(z)
    cof = [76.18009172947146; -86.50532032941677; 24.01409824083091; ...
           -1.231739572450155; 0.1208650973866179e-2; -0.5395239384953e-5];
    zp = z + (1:6)';
    ser = sum(cof ./ zp) + 1.000000000190015;
    tmp = z + 5.5 - (z + 0.5) * log(z + 5.5);
    f = -tmp + log(2.5066282746310005 * ser / z);
end
