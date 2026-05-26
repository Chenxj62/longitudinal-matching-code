function [total_M12, total_M22, total_M34, total_M44, k1_val] = M12opt(beamline_fwd, E0)

    me_c2   = 0.511e6;      % eV
    c_light = 299792458;
    num_elems = length(beamline_fwd);

    gamma_in_prof_fwd  = zeros(num_elems, 1);
    gamma_out_prof_fwd = zeros(num_elems, 1);
    family_G_dict = containers.Map();

    % =====================================================================
    % 1. 正向物理通道：计算能量剖面 & 锚定真实 k1
    % =====================================================================
    current_gamma = E0 / me_c2;
    quad_idx = 0;
    k1_val = [];

    for i = 1 : num_elems
        gamma_in_prof_fwd(i) = current_gamma;
        elem = beamline_fwd(i);

        if strcmpi(elem.type, 'cavity')
            if isfield(elem, 'num_cells') && isfield(elem, 'frequency')
                elem.length = elem.num_cells * 0.5 * c_light / (elem.frequency * 1e6);
                beamline_fwd(i).length = elem.length;
            end

            A_tmp   = 0; if isfield(elem, 'gradient'), A_tmp   = elem.gradient * 1e6; end
            phi_tmp = 0; if isfield(elem, 'phase'),    phi_tmp = elem.phase * pi/180; end

            current_gamma = current_gamma + (A_tmp * elem.length * cos(phi_tmp)) / me_c2;

        elseif strcmpi(elem.type, 'quadrupole')
            qname = '';
            if isfield(elem, 'name'), qname = elem.name; end

            if ~isempty(qname)
                quad_idx = quad_idx + 1;

                if family_G_dict.isKey(qname)
                    G_norm = family_G_dict(qname);
                    elem.k1 = G_norm / current_gamma;
                else
                    G_norm = elem.k1 * current_gamma;
                    family_G_dict(qname) = G_norm;
                end

                k1_val(quad_idx) = elem.k1;
                beamline_fwd(i) = elem;
            end
        end

        gamma_out_prof_fwd(i) = current_gamma;
    end

    % =====================================================================
    % 2. 数组翻转通道：调用外部翻转，对齐能量
    % =====================================================================
    beamline_rev = rev_beamline(beamline_fwd);

    gamma_in_prof_rev  = flip(gamma_in_prof_fwd);
    gamma_out_prof_rev = flip(gamma_out_prof_fwd);
    gamma_final        = gamma_out_prof_fwd(end);

    % =====================================================================
    % 3. 伴随矩阵逆向提取通道：四个矩阵元同步倒退 (R \ V_mat)
    %
    %    列1: x 平面，用 [0; -1] 提取 M12
    %    列2: x 平面，用 [1;  0] 提取 M22
    %    列3: y 平面，用 [0; -1] 提取 M34
    %    列4: y 平面，用 [1;  0] 提取 M44
    %
    % 对于 2x2 横向块 R = [a b; c d]，
    %   R^{-1}[0; -1] = [ b; -a ] / det(R)
    %   R^{-1}[1;  0] = [ d; -c ] / det(R)
    %
    % 所以在 cavity 中点：
    %   V_mid(1,1) * det = M12
    %   V_mid(1,2) * det = M22
    %   V_mid(3,3) * det = M34
    %   V_mid(3,4) * det = M44
    %
    % 当前 det(R_mid->final) = gamma_mid / gamma_final
    % =====================================================================
    V_acc = zeros(6, 4);

    % x-plane
    V_acc(2, 1) = -1;   % for M12
    V_acc(1, 2) =  1;   % for M22

    % y-plane
    V_acc(4, 3) = -1;   % for M34
    V_acc(3, 4) =  1;   % for M44

    total_M12 = 0;
    total_M22 = 0;
    total_M34 = 0;
    total_M44 = 0;

    for i = 1 : num_elems
        elem = beamline_rev(i);
        gamma_in  = gamma_in_prof_rev(i);
        gamma_out = gamma_out_prof_rev(i);

        if strcmpi(elem.type, 'cavity')
            [R_half1, R_half2] = getCavitySplitMatrixRS(elem, gamma_in, gamma_out, me_c2);

            % 退后半步到达 cavity 中点
            V_mid = R_half2 \ V_acc;

            % 计算 cavity 中点能量
      
          
   
                


            scale = 1/ gamma_final;

            % 恢复真实矩阵元
            M12_center = V_mid(1,1) * scale;
            M22_center = V_mid(1,2) * scale;
            M34_center = V_mid(3,3) * scale;
            M44_center = V_mid(3,4) * scale;

            % 如需“绝对值求和”，把下面四行改成 total = total + abs(...)
            total_M12 = total_M12 + M12_center;
            total_M22 = total_M22 + M22_center;
            total_M34 = total_M34 + M34_center;
            total_M44 = total_M44 + M44_center;

            % 继续退前半步到达 cavity 入口
            V_acc = R_half1 \ V_mid;

        else
            % 常规元件同步倒推
            R_elem = getElementMatrix(elem);
            V_acc = R_elem \ V_acc;
        end
    end

end


function [R_half1, R_half2] = getCavitySplitMatrixRS(elem, gamma_i, gamma_f, me_c2) %#ok<INUSD>
    L = elem.length;

    R_half1 = eye(6);
    R_half2 = eye(6);

    if L > 0
        gamma_prime = (gamma_f - gamma_i) / L;

        if abs(gamma_prime) > 1e-10
            Omega = 1 / sqrt(8);
            gamma_mid = gamma_i + gamma_prime * (L / 2);

            % -------------------------
            % R_half1: 入口 Kick + 前半段 Body
            % -------------------------
            R_in = [1, 0; ...
                   -gamma_prime/(2*gamma_i), 1];

            ps1 = Omega * log(gamma_mid / gamma_i);

            R_body1 = [cos(ps1), (gamma_i/(gamma_prime*Omega))*sin(ps1); ...
                      -(gamma_prime*Omega/gamma_mid)*sin(ps1), (gamma_i/gamma_mid)*cos(ps1)];

            R_trans1 = R_body1 * R_in;

            R_half1(1:2,1:2) = R_trans1;
            R_half1(3:4,3:4) = R_trans1;
            R_half1(6,6) = gamma_i / gamma_mid;

            % -------------------------
            % R_half2: 后半段 Body + 出口 Kick
            % -------------------------
            ps2 = Omega * log(gamma_f / gamma_mid);

            R_body2 = [cos(ps2), (gamma_mid/(gamma_prime*Omega))*sin(ps2); ...
                      -(gamma_prime*Omega/gamma_f)*sin(ps2), (gamma_mid/gamma_f)*cos(ps2)];

            R_out = [1, 0; ...
                     gamma_prime/(2*gamma_f), 1];

            R_trans2 = R_out * R_body2;

            R_half2(1:2,1:2) = R_trans2;
            R_half2(3:4,3:4) = R_trans2;
            R_half2(6,6) = gamma_mid / gamma_f;

        else
            % 近似无加速时，退化为两段半漂移
            R_half1(1,2) = L/2;
            R_half1(3,4) = L/2;

            R_half2(1,2) = L/2;
            R_half2(3,4) = L/2;
        end
    else
        % 零长度保护
        R_half1 = eye(6);
        R_half2 = eye(6);
    end
end