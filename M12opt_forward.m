function [total_M12, total_M34, k1_val] = M12opt_forward(beamline_fwd, E0)
% M12opt_forward  Forward-pass benchmark for M12opt (adjoint version).
%
% Computes identical total_M12, total_M34, k1_val without reversing the
% beamline or using the adjoint method.
%
% Key identity:
%   In M12opt, the adjoint variable at cavity c's midpoint satisfies
%       V_mid(1,1) = R_mid_to_end(1,2)
%       V_mid(3,2) = R_mid_to_end(3,4)
%   where R_mid_to_end = transfer matrix from that midpoint to the exit.
%
% This benchmark builds R_mid_to_end by accumulating a "right product"
% R_right backward from the beamline end:
%
%     R_right = R_N * R_{N-1} * ... * R_{i+1}   (product of all elements after i)
%
% Recurrence (i going from N down to 1):
%     R_right_new = R_right_old * R_i
%
% For cavity i:
%     R_mid_to_end = R_right * R_half2          (NO inversion needed)
%
% This avoids the numerical errors of R_total * inv(R_to_mid).

    me_c2   = 0.511e6;      % eV
    c_light = 299792458;    % m/s
    num_elems = length(beamline_fwd);

    % =====================================================================
    % 1. Forward energy pass: gamma profile + k1 anchoring
    %    (identical to M12opt)
    % =====================================================================
    current_gamma  = E0 / me_c2;
    gamma_in_prof  = zeros(num_elems, 1);
    gamma_out_prof = zeros(num_elems, 1);
    family_G_dict  = containers.Map();
    quad_idx = 0;
    k1_val   = [];

    for i = 1 : num_elems
        gamma_in_prof(i) = current_gamma;
        elem = beamline_fwd(i);

        if strcmpi(elem.type, 'cavity')
            if isfield(elem, 'num_cells') && isfield(elem, 'frequency')
                elem.length = elem.num_cells * 0.5 * c_light / (elem.frequency * 1e6);
                beamline_fwd(i).length = elem.length;
            end
            A_tmp   = 0; if isfield(elem, 'gradient'), A_tmp   = elem.gradient * 1e6; end
            phi_tmp = 0; if isfield(elem, 'phase'),    phi_tmp = elem.phase * pi/180;  end
            current_gamma = current_gamma + (A_tmp * elem.length * cos(phi_tmp)) / me_c2;

        elseif strcmpi(elem.type, 'quadrupole')
            qname = ''; if isfield(elem, 'name'), qname = elem.name; end
            quad_idx = quad_idx + 1;
            if ~isempty(qname)
                if family_G_dict.isKey(qname)
                    elem.k1 = family_G_dict(qname) / current_gamma;
                else
                    family_G_dict(qname) = elem.k1 * current_gamma;
                end
                k1_val(quad_idx) = elem.k1;
                beamline_fwd(i)  = elem;
            end
        end

        gamma_out_prof(i) = current_gamma;
    end

    gamma_final = gamma_out_prof(end);

    % =====================================================================
    % 2. Forward pass: compute and cache every element's matrix.
    %    For cavities, cache R_half1, R_half2, and gamma_mid.
    % =====================================================================
    elem_R      = cell(num_elems, 1);   % for non-cavity elements
    cav_R_half1 = cell(num_elems, 1);   % for cavities
    cav_R_half2 = cell(num_elems, 1);
    is_cav      = false(num_elems, 1);

    for i = 1 : num_elems
        elem      = beamline_fwd(i);
        gamma_in  = gamma_in_prof(i);
        gamma_out = gamma_out_prof(i);

        if strcmpi(elem.type, 'cavity')
            [R_half1, R_half2] = getCavitySplitMatrixRS(elem, gamma_in, gamma_out, me_c2);
            is_cav(i)      = true;
            cav_R_half1{i} = R_half1;
            cav_R_half2{i} = R_half2;
        else
            elem_R{i} = getElementMatrix(elem);
        end
    end

    % =====================================================================
    % 3. Backward accumulation pass.
    %
    %    R_right starts as I and grows rightward as i decreases:
    %      R_right = R_N * R_{N-1} * ... * R_{i+1}
    %
    %    For cavity at position i:
    %      R_mid_to_end = R_right * R_half2
    %                   = (product of all elements after i) * R_half2
    %                   = transfer matrix from cavity midpoint to beamline end
    %
    %    Update: R_right = R_right * R_half2 * R_half1  (include full cavity)
    %    For non-cavity: R_right = R_right * R_elem
    % =====================================================================
    R_right   = eye(6);
    total_M12 = 0;
    total_M34 = 0;

    for i = num_elems : -1 : 1
        if is_cav(i)
            R_half1 = cav_R_half1{i};
            R_half2 = cav_R_half2{i};

            R_mid_to_end = R_right * R_half2;   % no inversion needed

            % No gamma weight needed: in the adjoint method the explicit
            % (gamma_mid/gamma_final) compensates the 1/det factor from
            % inverting the non-symplectic matrix (det_x = gamma_mid/gamma_final).
            % Here we read R(1,2) directly — adiabatic damping is already
            % encoded in the gamma-dependent RS matrix elements.
            total_M12 = total_M12 + R_mid_to_end(1, 2);
            total_M34 = total_M34 + R_mid_to_end(3, 4);

            R_right = R_right * R_half2 * R_half1;
        else
            R_right = R_right * elem_R{i};
        end
    end

end


% -------------------------------------------------------------------------
% Local copy of getCavitySplitMatrixRS (same as in M12opt.m)
% -------------------------------------------------------------------------
function [R_half1, R_half2] = getCavitySplitMatrixRS(elem, gamma_i, gamma_f, me_c2)
    L = elem.length;
    gamma_prime = (gamma_f - gamma_i) / L;
    R_half1 = eye(6); R_half2 = eye(6);

    if L > 0 && abs(gamma_prime) > 1e-10
        Omega     = 1 / sqrt(8);
        gamma_mid = gamma_i + gamma_prime * (L / 2);

        R_in    = [1, 0; -gamma_prime/(2*gamma_i), 1];
        ps1     = Omega * log(gamma_mid / gamma_i);
        R_body1 = [cos(ps1),  (gamma_i/(gamma_prime*Omega))*sin(ps1); ...
                  -(gamma_prime*Omega/gamma_mid)*sin(ps1), (gamma_i/gamma_mid)*cos(ps1)];
        R_trans1 = R_body1 * R_in;
        R_half1(1:2, 1:2) = R_trans1; R_half1(3:4, 3:4) = R_trans1;
        R_half1(6,6) = gamma_i / gamma_mid;

        ps2     = Omega * log(gamma_f / gamma_mid);
        R_body2 = [cos(ps2),  (gamma_mid/(gamma_prime*Omega))*sin(ps2); ...
                  -(gamma_prime*Omega/gamma_f)*sin(ps2), (gamma_mid/gamma_f)*cos(ps2)];
        R_out   = [1, 0; gamma_prime/(2*gamma_f), 1];
        R_trans2 = R_out * R_body2;
        R_half2(1:2, 1:2) = R_trans2; R_half2(3:4, 3:4) = R_trans2;
        R_half2(6,6) = gamma_mid / gamma_f;
    else
        R_half1(1,2) = L/2; R_half1(3,4) = L/2;
        R_half2(1,2) = L/2; R_half2(3,4) = L/2;
    end
end
