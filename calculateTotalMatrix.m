function [R_total, gamma_out] = calculateTotalMatrix(beamline, gamma_in)
% calculateTotalMatrix  Compute the cumulative 6x6 first-order transport matrix
%
% Usage:
%   R = calculateTotalMatrix(beamline)
%       Original behaviour: no energy tracking, gamma not required.
%
%   [R, gamma_out] = calculateTotalMatrix(beamline, gamma_in)
%       Track Lorentz gamma through the beamline (needed for cavities).
%       gamma_out = gamma at the beamline exit.

    R_total   = eye(6);
    gamma_out = [];

    use_gamma = (nargin >= 2) && ~isempty(gamma_in);
    if use_gamma
        gamma_out = gamma_in;
    end

    for i = 1:length(beamline)
        if use_gamma
            [R_elem, gamma_out] = getElementMatrix(beamline(i), gamma_out);
        else
            R_elem = getElementMatrix(beamline(i));
        end
        R_total = R_elem * R_total;
    end
end

