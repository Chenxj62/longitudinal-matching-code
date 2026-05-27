function ds = local_s_in_dipole(s, dip)
% local_s_in_dipole  Distance from dipole entrance [cm]
%
% Returns how far inside the current dipole the position s is.
% Returns 0 if s is not inside any dipole.
%
% Inputs:
%   s   - position [m]
%   dip - dipole table struct (.N, .s1, .s2, .rho) all [m]

    ds = 0;
    for m = 1:dip.N
        if s >= dip.s1(m) && s <= dip.s2(m)
            ds = s - dip.s1(m);
        end
    end
end
