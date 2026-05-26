function [R, gamma_out] = getElementMatrix(element, gamma)
    % 计算元件的传输矩阵
    % 输入:
    %   element - 元件结构体
    %   gamma - 相对论因子 (可选，仅cavity需要)
    % 输出:
    %   R - 6x6传输矩阵
    %   gamma_out - 出口处相对论因子 (可选，非cavity元件原样返回)
    
    R = eye(6);
    L = element.length;
    
    % 默认出口gamma等于入口gamma
    if nargin >= 2 && ~isempty(gamma)
        gamma_out = gamma;
    else
        gamma_out = [];
    end
    
    % 物理常数
    m0 = 0.511;  % 电子静止质量 (MeV)
    
    switch lower(element.type)
        case 'edge'
            h = element.h;
            angle = element.angle;
            R(2,1) = h*tan(angle);
            R(4,3) = -R(2,1);
            
        case 'drift'
            R(1,2) = L;
            R(3,4) = L;
            
        case 'dipole'
            h = element.h;
            if abs(h) < 1e-10
                R(1,2) = L;
                R(3,4) = L;
            else
                theta = L * h;
                cos_t = cos(theta);
                sin_t = sin(theta);
                R(1,1:2) = [cos_t, sin_t/h];
                R(3,4) = L;
                R(1,6) = (1 - cos_t)/h;
                R(2,1:2) = [-h*sin_t, cos_t];
                R(2,6) = sin_t;
                R(5,1) = -sin_t;
                R(5,2) = -(1 - cos_t)/h;
                R(5,6) = -(theta - sin_t)/h;
            end
            
        case {'quadrupole', 'quad'}
            k_actual = element.k1;
            if abs(k_actual) < 1e-10
                R(1,2) = L;
                R(3,4) = L;
            else
                if k_actual > 0
                    % 水平聚焦，垂直散焦
                    kl = sqrt(k_actual) * L;
                    R(1,1) = cos(kl);
                    R(1,2) = sin(kl)/sqrt(k_actual);
                    R(2,1) = -sqrt(k_actual)*sin(kl);
                    R(2,2) = cos(kl);
                    
                    R(3,3) = cosh(kl);
                    R(3,4) = sinh(kl)/sqrt(k_actual);
                    R(4,3) = sqrt(k_actual)*sinh(kl);
                    R(4,4) = cosh(kl);
                elseif k_actual < 0
                    % 水平散焦，垂直聚焦
                    kl = sqrt(-k_actual) * L;
                    R(1,1) = cosh(kl);
                    R(1,2) = sinh(kl)/sqrt(-k_actual);
                    R(2,1) = sqrt(-k_actual)*sinh(kl);
                    R(2,2) = cosh(kl);
                    
                    R(3,3) = cos(kl);
                    R(3,4) = sin(kl)/sqrt(-k_actual);
                    R(4,3) = -sqrt(-k_actual)*sin(kl);
                    R(4,4) = cos(kl);
                end
            end
            
        case {'sextupole', 'sext'}
            R(1,2) = L;
            R(3,4) = L;
            
        case 'marker'
            % Marker不改变矩阵，直接返回 eye(6)
            
        case {'cavity', 'rfcavity', 'rf'}
            % ==========================================
            % 重构：基于严格 Rosenzweig-Serafini 理论
            % ==========================================
            if nargin < 2 || isempty(gamma)
                error('计算cavity传输矩阵需要提供gamma（相对论因子）参数');
            end
            
            gradient  = 0; if isfield(element, 'gradient'), gradient = element.gradient; end
            phase_deg = 0; if isfield(element, 'phase'),    phase_deg = element.phase; end
            Delta_phi = phase_deg * pi/180;
            
            % 计算能量增益 (MV)
            if isfield(element, 'voltage')
                energy_gain = element.voltage * cos(Delta_phi);
            else
                energy_gain = gradient * L * cos(Delta_phi);
            end
            
            gamma_i = gamma;
            gamma_f = gamma_i + energy_gain / m0;
            
            if gamma_f <= 0
                warning('Cavity: 出口能量非正值，设置为入口值');
                gamma_f = gamma_i;
            end
            gamma_out = gamma_f;
            gamma_prime = (gamma_f - gamma_i) / L;
            
            % --- RS 横向传输矩阵构建 ---
            if L > 0 && abs(gamma_prime) > 1e-10
                Omega = 1 / sqrt(8);  % pi-mode standing wave 聚焦因子
                phi_phase = Omega * log(gamma_f / gamma_i); % Betatron 相移
                
                % 1. 入口边缘场 (Entrance Kick)
                R_in = [1, 0; -gamma_prime/(2*gamma_i), 1];
                
                % 2. 腔体主体 (Cavity Body)
                R_body = [cos(phi_phase),                      (gamma_i/(gamma_prime*Omega))*sin(phi_phase); ...
                         -(gamma_prime*Omega/gamma_f)*sin(phi_phase), (gamma_i/gamma_f)*cos(phi_phase)];
                     
                % 3. 出口边缘场 (Exit Kick)
                R_out = [1, 0; gamma_prime/(2*gamma_f), 1];
                
                % 组合总横向矩阵
                R_trans = R_out * R_body * R_in;
                
                % 赋值给 x 和 y 平面
                R(1:2, 1:2) = R_trans;
                R(3:4, 3:4) = R_trans; 
                R(4,4)=1*R(4,4);
            else
                % 无加速时退化为漂移段
                R(1,2) = L;
                R(3,4) = L;
            end
            
            % --- 纵向传输矩阵 ---
            if isfield(element, 'frequency')
                krf = 2 * pi * element.frequency / 299.792458;
                R(6,5) = (1 - gamma_i / gamma_f) * krf * tan(Delta_phi);
            end
            R(5,5) = 1;
            R(6,6) = gamma_i / gamma_f;
            
        case 'tgu' 
            gam = element.gamma; 
            alp = element.alp;   
            alpc = element.alpc;  
            K0 = element.k0;
            ku = element.ku;
            if abs(alp - alpc) < 1e-10
                R(1,2) = L;
                R(3,4) = L;
            else
                eta_c = 1 / (alp - alpc);
                eta_alpha = (2 + K0^2) / (alp * K0^2);
                
                kx = sqrt((alp * K0^2) / (2 * gam^2) * (alp - alpc));
                ky = sqrt(K0^2 / (2 * gam^2) * (ku^2 + alp^2 + alp * alpc));
                
                if (alp - alpc) < 0
                    warning('TGU: alp - alpc < 0, kx becomes imaginary. Check parameters.');
                end
                
                cos_kxz = cos(kx * L); sin_kxz = sin(kx * L);
                cos_kyz = cos(ky * L); sin_kyz = sin(ky * L);
                
                R(1,1) = cos_kxz; R(1,2) = (1/kx) * sin_kxz; R(1,6) = eta_c * (1 - cos_kxz);
                R(2,1) = -kx * sin_kxz; R(2,2) = cos_kxz; R(2,6) = eta_c * kx * sin_kxz;
                R(3,3) = cos_kyz; R(3,4) = (1/ky) * sin_kyz;
                R(4,3) = -ky * sin_kyz; R(4,4) = cos_kyz;
                R(5,1) = -eta_c * kx * sin_kxz; R(5,2) = -eta_c * (1 - cos_kxz);
                R(5,6) = kx^2 * eta_c * eta_alpha * L - eta_c^2 * kx * (kx * L - sin_kxz);
            end
            
        case '1stttay'
            R(1,1) = element.m11; R(2,2) = element.m22;
            R(2,1) = element.m21; R(1,2) = ( element.m11 * element.m22-1) / element.m21;
            R(3,3) = element.m33; R(4,4) = element.m44;
            R(4,3) = element.m43; R(3,4) = ( element.m33 * element.m44-1) / element.m43;
            
        otherwise
            R(1,2) = L;
            R(3,4) = L;
    end
end