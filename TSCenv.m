function [results] = TSCenv(beamline, E0, beamfile1)
% TSCenv - 空间电荷包络计算（加入 drift 中的速度压缩对 Ksc 的影响）
%
% 输入：
%   beamline  : 束线元件结构体数组
%   E0        : 参考能量 [eV]
%   beamfile1 : 切片参数文件 (由 read_z_distribution 生成)
%               第一行 col16 = Q_total [C]

% --- 读取切片文件获取全局参数 ---
try
    slice_data_matrix_pre = load(beamfile1);
catch
    error('无法读取文件: %s', beamfile1);
end

Q_total = slice_data_matrix_pre(1, 16)*1/1.1;   % bunch charge [C]

calc_k1_flag = 1;
me_c2   = 0.511e6;
c_light = 2.998e8;
ds_kick = 0.02;

% 速度压缩 R56 = -L/(beta^2 * gamma^2)，见 driftR56_velocity()

% ========================================================================
% 1. 提前获取参考能量 E_ref 用于全局预处理，确保全线归一化标定统一
% ========================================================================
z_cent_all = slice_data_matrix_pre(2:end-1, 14);
delta_all  = slice_data_matrix_pre(2:end-1, 15);
[~, min_z_idx] = min(abs(z_cent_all));
E_ref_global = E0 * (1 + delta_all(min_z_idx));

% ========================================================================
% 2. 全局预处理：绑定同名四极铁的物理梯度 G_norm，统一腔长
% ========================================================================
temp_E = E_ref_global;
quad_name_map = containers.Map();

for i = 1:length(beamline)
    gamma_entrance = temp_E / me_c2;
    elem = beamline(i);

    if strcmpi(elem.type, 'cavity')
        if isfield(elem, 'num_cells') && isfield(elem, 'frequency')
            beamline(i).length = elem.num_cells * 0.5 * c_light / (elem.frequency * 1e6);
        end
        A_tmp   = 0; if isfield(elem, 'gradient'), A_tmp   = elem.gradient * 1e6; end
        phi_tmp = 0; if isfield(elem, 'phase'),    phi_tmp = elem.phase * pi/180; end
        temp_E  = temp_E + A_tmp * beamline(i).length * cos(phi_tmp);

    elseif strcmpi(elem.type, 'quadrupole')
        qname = '';
        if isfield(elem, 'name'), qname = elem.name; end
        G_norm_here = elem.k1 * gamma_entrance;

        if ~isempty(qname) && quad_name_map.isKey(qname)
            bound_G_norm = quad_name_map(qname);
            beamline(i).k1 = bound_G_norm / gamma_entrance;
        else
            if ~isempty(qname)
                quad_name_map(qname) = G_norm_here;
            end
        end
    end
end

% ================= 基于切片的多包络计算模式 =================
slice_data_matrix = slice_data_matrix_pre;

total_particles = slice_data_matrix(1, 1);

slice_info      = slice_data_matrix(2:end-1, :);
num_slice       = size(slice_info, 1);

slice_particle_counts = slice_info(:, 1);
betax_slices          = slice_info(:, 2);
alphax_slices         = slice_info(:, 3);
emitx_geo_slices      = slice_info(:, 5);
betay_slices          = slice_info(:, 6);
alphay_slices         = slice_info(:, 7);
emity_geo_slices      = slice_info(:, 9);
x_centroid_slices     = slice_info(:, 10);
px_centroid_slices    = slice_info(:, 11);
y_centroid_slices     = slice_info(:, 12);
py_centroid_slices    = slice_info(:, 13);
z_centroid_slices     = slice_info(:, 14);
delta_centroid_slices = slice_info(:, 15);
x_pz_slices           = slice_info(:, 16);
px_pz_slices          = slice_info(:, 17);

rc = 2.8179e-15;
e  = 1.602e-19;

z_val0 = z_centroid_slices(2)-z_centroid_slices(1);

if calc_k1_flag
    k1_values = [];

end

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep', ds_kick, 'Stats', 'off');

valid_slices = slice_particle_counts > 0 & emitx_geo_slices > 1e-15 & ...
               emity_geo_slices > 1e-15 & betax_slices > 1e-10 & betay_slices > 1e-10;
valid_slice_indices = find(valid_slices);

% --- 确认参考切片能量 ---
z_valid = z_centroid_slices(valid_slice_indices);
[~, ref_local_idx] = min(abs(z_valid));
ref_slice_idx = valid_slice_indices(ref_local_idx);
E_ref         = E0 * (1 + delta_centroid_slices(ref_slice_idx));
gamma_ref     = E_ref / me_c2;

% --- 用全束首尾切片估计线性 chirp: delta ≈ h_delta_global * z ---
if numel(valid_slice_indices) >= 2
    [z_sorted, ord] = sort(z_centroid_slices(valid_slice_indices));
    d_sorted = delta_centroid_slices(valid_slice_indices);
    d_sorted = d_sorted(ord);

    dz_full = z_sorted(end) - z_sorted(1);
    if abs(dz_full) > 1e-15
        h_delta_global = (d_sorted(end) - d_sorted(1)) / dz_full;
    else
        h_delta_global = 0;
    end
else
    h_delta_global = 0;
end

% --- 预计算非线性项 (可选) ---
ilist = [1, 1, 2, 2, 3, 3, 4, 4];
jlist = [2, 1, 1, 2, 4, 3, 3, 4];
klist = [6, 6, 6, 6, 6, 6, 6, 6];
[s_transport, T_evol] = nonlinear_calculator(beamline, ilist, jlist, klist, 0.001, E_ref);
T126_evol = T_evol(:,1);  T116_evol = T_evol(:,2);
T216_evol = T_evol(:,3);  T226_evol = T_evol(:,4);
T346_evol = T_evol(:,5);  T336_evol = T_evol(:,6);
T436_evol = T_evol(:,7);  T446_evol = T_evol(:,8);

num_elements     = length(beamline);
total_length     = sum([beamline.length]);
estimated_points = ceil(total_length / ds_kick) + num_elements * 5;

all_slice_envelope(num_slice) = struct( ...
    's', [], 'sigx', [], 'sigpx', [], 'sigy', [], 'sigpy', [], ...
    'E', [], 'epsn2x', [], 'epsn2y', [], 'particle_count', 0, 'valid', false);

all_slice_centroid(num_slice) = struct( ...
    's_exit', [], 'x_exit', [], 'px_exit', [], 'y_exit', [], 'py_exit', [], ...
    'particle_count', 0, 'valid', false);

s_envelope_common = [];
element_types = zeros(num_elements, 1);

for i = 1:num_elements
    switch lower(beamline(i).type)
        case 'drift'
            element_types(i) = 1;
        case 'quadrupole'
            element_types(i) = 2;
        case 'cavity'
            element_types(i) = 3;
        case {'sextupole','marker'}
            element_types(i) = 0;
        otherwise
            element_types(i) = 4;
    end
end

first_cavity_idx = find(element_types == 3, 1, 'first');
if isempty(first_cavity_idx)
    first_cavity_idx = inf;
end

% --- 预计算参考能量演化 (用于 Delta_eff 计算) ---
E_ref_at_element = zeros(num_elements, 1);
E_running_ref = E_ref;
for i = 1:num_elements
    if element_types(i) == 3
        A_r   = 0; if isfield(beamline(i), 'gradient'), A_r = beamline(i).gradient * 1e6; end
        phi_r = 0; if isfield(beamline(i), 'phase'),    phi_r = beamline(i).phase * pi/180; end
        E_running_ref = E_running_ref + A_r * beamline(i).length * cos(phi_r);
    end
    E_ref_at_element(i) = E_running_ref;
end

% ================= 逐切片极速计算核心 =================
for idx = 1:length(valid_slice_indices)
    slice_idx = valid_slice_indices(idx);

    betax_slice      = betax_slices(slice_idx);   alphax_slice = alphax_slices(slice_idx);
    betay_slice      = betay_slices(slice_idx);   alphay_slice = alphay_slices(slice_idx);
    emitx_geo_slice  = emitx_geo_slices(slice_idx);
    emity_geo_slice  = emity_geo_slices(slice_idx);
    slice_particles  = slice_particle_counts(slice_idx);

    x0_slice  = x_centroid_slices(slice_idx);  px0_slice  = px_centroid_slices(slice_idx);
    y0_slice  = y_centroid_slices(slice_idx);  py0_slice  = py_centroid_slices(slice_idx);
    delta0_slice = delta_centroid_slices(slice_idx);
    z0_slice     = z_centroid_slices(slice_idx);

    sig0x_slice  = sqrt(betax_slice * emitx_geo_slice);
    sigp0x_slice = -alphax_slice * sqrt(emitx_geo_slice / betax_slice);
    sig0y_slice  = sqrt(betay_slice * emity_geo_slice);
    sigp0y_slice = -alphay_slice * sqrt(emity_geo_slice / betay_slice);

    x_pz_slice  = x_pz_slices(slice_idx);
    px_pz_slice = px_pz_slices(slice_idx);

    if ~all(isfinite([sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice])) || ...
       sig0x_slice <= 0 || sig0y_slice <= 0
        continue;
    end

    Q_slice      = Q_total * slice_particles / total_particles;
    N_slice_real = abs(Q_slice) / e;

    E_slice_init     = E0 * (1 + delta0_slice);
    gamma_slice_init = E_slice_init / me_c2;
    epsn2x0 = (emitx_geo_slice * gamma_slice_init)^2;
    epsn2y0 = (emity_geo_slice * gamma_slice_init)^2;

    current_y_slice = [sig0x_slice, sigp0x_slice, sig0y_slice, sigp0y_slice, E_slice_init, epsn2x0, epsn2y0]';
    current_s_slice = 0;

    % 当前 slice 的 z_val：只在第一个 cavity 前的 drift 中演化
    z_val_current = z_val0;
    freeze_zval_after_first_cavity = false;

    % ---------- 质心矩阵跟踪 ----------
    num_exit_points = length(s_transport);
    xc = zeros(num_exit_points, 1); pxc = zeros(num_exit_points, 1);
    yc = zeros(num_exit_points, 1); pyc = zeros(num_exit_points, 1);
    xc(1) = x0_slice; pxc(1) = px0_slice; yc(1) = y0_slice; pyc(1) = py0_slice;

    for ii = 2:num_exit_points
        R = eye(6);
        temp_gam = gamma_ref * (1 + delta0_slice);
        temp_E_slice = temp_gam * me_c2;
        for j = 1:ii-1
            elem = beamline(j);
            if strcmpi(elem.type, 'cavity')
                R_elem  = getElementMatrix(elem, temp_gam);
                A_c     = 0; if isfield(elem, 'gradient'), A_c = elem.gradient * 1e6; end
                phi_c   = 0; if isfield(elem, 'phase'),    phi_c = elem.phase * pi/180; end
                k_RF_c  = 2*pi*elem.frequency*1e6 / c_light;
                phi_eff_c = phi_c + k_RF_c * z0_slice;
                temp_E_slice = temp_E_slice + A_c * elem.length * cos(phi_eff_c);
                temp_gam = temp_E_slice / me_c2;
            elseif strcmpi(elem.type, 'quadrupole') || strcmpi(elem.type, 'quad')
                elem_chr = elem;
                delta_eff_c = (temp_E_slice - E_ref_at_element(j)) / E_ref_at_element(j);
                elem_chr.k1 = elem.k1 / (1 + delta_eff_c);
                R_elem = getElementMatrix(elem_chr);
            else
                R_elem = getElementMatrix(elem);
            end
            R = R_elem * R;
        end
        xc(ii)  =R(1,1)*x0_slice + R(1,2)*px0_slice+ T116_evol(ii) * x_pz_slice  + T126_evol(ii) * px_pz_slice;
        pxc(ii) =R(2,1)*x0_slice + R(2,2)*px0_slice + T216_evol(ii) * x_pz_slice  + T226_evol(ii) * px_pz_slice;
        yc(ii)  =R(3,3)*y0_slice + R(3,4)*py0_slice;
        pyc(ii) =R(4,3)*y0_slice + R(4,4)*py0_slice;
    end

    % ---------- 包络跟踪 ----------
    s_env      = zeros(estimated_points, 1);
    sigx_env   = zeros(estimated_points, 1); sigpx_env  = zeros(estimated_points, 1);
    sigy_env   = zeros(estimated_points, 1); sigpy_env  = zeros(estimated_points, 1);
    E_env      = zeros(estimated_points, 1);
    epsn2x_env = zeros(estimated_points, 1); epsn2y_env = zeros(estimated_points, 1);

    cur_idx = 0;
    ode_failed = false;

    for iel = 1:num_elements
        element = beamline(iel);
        if element.length == 0 || element_types(iel) == 0
            continue;
        end

        s_start = current_s_slice;
        s_end   = current_s_slice + element.length;

        if any(~isfinite(current_y_slice)) || current_y_slice(1) <= 0 || current_y_slice(3) <= 0
            ode_failed = true;
            break;
        end

        % 默认：该元件入口 z_val 固定
        z_val_elem_in = z_val_current;
        evolve_zval_flag = false;
        kz_drift = 0;

        % 只在第一个 cavity 之前的 drift 中做速度压缩
        if (iel < first_cavity_idx) && (element_types(iel) == 1) && ~freeze_zval_after_first_cavity
            gamma_here = max(current_y_slice(5) / me_c2, 1.0000001);
            R56v_total = driftR56_velocity(element.length, gamma_here);
            kz_drift   = h_delta_global * (R56v_total / max(element.length, 1e-15));
            evolve_zval_flag = true;
        end

        if current_y_slice(5) >= 100e6
            % --- 高能解析路径 (忽略空电) ---
            gamma_in     = current_y_slice(5) / me_c2;
            epsnx2_in    = current_y_slice(6);
            epsny2_in    = current_y_slice(7);
            emitx_geo_in = sqrt(max(epsnx2_in, 0)) / gamma_in;
            emity_geo_in = sqrt(max(epsny2_in, 0)) / gamma_in;

            sigx_in  = current_y_slice(1);  sigpx_in = current_y_slice(2);
            sigy_in  = current_y_slice(3);  sigpy_in = current_y_slice(4);

            Sigma_x_in = [sigx_in^2, sigx_in*sigpx_in; ...
                          sigx_in*sigpx_in, (emitx_geo_in^2/sigx_in^2) + sigpx_in^2];
            Sigma_y_in = [sigy_in^2, sigy_in*sigpy_in; ...
                          sigy_in*sigpy_in, (emity_geo_in^2/sigy_in^2) + sigpy_in^2];

            elem_mat = element;
            if element_types(iel) == 2
                k1_nominal = element.k1;
                delta_eff  = (current_y_slice(5) - E_ref_at_element(iel)) / E_ref_at_element(iel);
                elem_mat.k1 = k1_nominal / (1 + delta_eff);
                if calc_k1_flag && slice_idx == valid_slice_indices(1)
                    k1_values = [k1_values, k1_nominal]; %#ok<AGROW>
                end
                R_elem = getElementMatrix(elem_mat);
                E_out  = current_y_slice(5);

            elseif element_types(iel) == 3
                 A = 0; if isfield(element, 'gradient'), A = element.gradient * 1e6; end
                phi = 0; if isfield(element, 'phase'), phi = element.phase * pi/180; end
                k_RF = 2*pi*element.frequency*1e6 / c_light;
                phi_eff = phi + k_RF * z0_slice;
                elem_mat.phase=phi_eff;
                R_elem  = getElementMatrix(elem_mat, gamma_in);
               
                E_out = current_y_slice(5) + A * element.length * cos(phi_eff);

            else
                R_elem = getElementMatrix(elem_mat);
                E_out  = current_y_slice(5);
            end

            Rx = R_elem(1:2, 1:2);
            Ry = R_elem(3:4, 3:4);

            Sigma_x_out = Rx * Sigma_x_in * Rx';
            Sigma_y_out = Ry * Sigma_y_in * Ry';

            sigx_out  = sqrt(max(Sigma_x_out(1,1), 1e-24));
            sigpx_out = Sigma_x_out(1,2) / sigx_out;
            sigy_out  = sqrt(max(Sigma_y_out(1,1), 1e-24));
            sigpy_out = Sigma_y_out(1,2) / sigy_out;

            z_seg = [s_start; s_end];
            y_seg = [current_y_slice'; sigx_out, sigpx_out, sigy_out, sigpy_out, E_out, epsnx2_in, epsny2_in];
            current_y_slice = y_seg(end,:)';

        else
            % --- 低能 ODE 路径 ---
            switch element_types(iel)
                case 1
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced( ...
                        z, y, 0, 0, 0, N_slice_real, z_val_elem_in, rc, me_c2, ...
                        s_start, kz_drift, evolve_zval_flag), ...
                        [s_start, s_end], current_y_slice, options);

                case 2
                    k1_nominal  = element.k1;
                    delta_eff    = (current_y_slice(5) - E_ref_at_element(iel)) / E_ref_at_element(iel);
                    k1_effective = k1_nominal / (1 + delta_eff);
                    if calc_k1_flag && slice_idx == valid_slice_indices(1)
                        k1_values = [k1_values, k1_nominal]; %#ok<AGROW>
                    end
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced( ...
                        z, y, k1_effective, 0, 0, N_slice_real, z_val_elem_in, rc, me_c2, ...
                        s_start, 0, false), ...
                        [s_start, s_end], current_y_slice, options);

                case 3
                    A = 0; if isfield(element, 'gradient'), A = element.gradient * 1e6; end
                    phi = 0; if isfield(element, 'phase'), phi = element.phase * pi/180; end
                    k_RF = 2*pi*element.frequency*1e6 / c_light;
                    phi_eff = phi + k_RF * z0_slice;

                    % 边缘场入口 Kick
                    dg_entrance = A * cos(phi) / current_y_slice(5);
                    M21_ent = -dg_entrance / 2;
                    current_y_slice(2) = current_y_slice(2) + M21_ent * current_y_slice(1);
                    current_y_slice(4) = current_y_slice(4) + M21_ent * current_y_slice(3);

                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced( ...
                        z, y, 0, A, phi_eff, N_slice_real, z_val_elem_in, rc, me_c2, ...
                        s_start, 0, false), ...
                        [s_start, s_end], current_y_slice, options);

                    if isempty(z_seg) || any(any(~isfinite(y_seg)))
                        ode_failed = true;
                        break;
                    end

                    % 边缘场出口 Kick
                    current_y_slice = y_seg(end,:)';
                    dg_exit  = A * cos(phi) / current_y_slice(5);
                    M21_exit = dg_exit / 2;
                    current_y_slice(2) = current_y_slice(2) + M21_exit * current_y_slice(1);
                    current_y_slice(4) = current_y_slice(4) + M21_exit * current_y_slice(3);

                otherwise
                    [z_seg, y_seg] = ode45(@(z,y) envelopeODE_sliced( ...
                        z, y, 0, 0, 0, N_slice_real, z_val_elem_in, rc, me_c2, ...
                        s_start, 0, false), ...
                        [s_start, s_end], current_y_slice, options);
            end

            if isempty(z_seg) || any(any(~isfinite(y_seg)))
                ode_failed = true;
                break;
            end

            if element_types(iel) ~= 3
                current_y_slice = y_seg(end,:)';
            end
        end

        % drift 结束后更新 z_val_current；进入第一个 cavity 后冻结
        if evolve_zval_flag
            z_val_current = z_val_elem_in * localZScale(element.length, kz_drift);
        end

        if iel == first_cavity_idx
            freeze_zval_after_first_cavity = true;
        end

        % [极速拼接方案]：切掉第一个点，保证 s 轴严格递增
        n_new = length(z_seg);
        if cur_idx == 0
            idx_range = 1:n_new;
        else
            idx_range = cur_idx+1 : cur_idx+n_new-1;
            z_seg = z_seg(2:end);
            y_seg = y_seg(2:end,:);
        end

        s_env(idx_range)      = z_seg;
        sigx_env(idx_range)   = y_seg(:,1);
        sigpx_env(idx_range)  = y_seg(:,2);
        sigy_env(idx_range)   = y_seg(:,3);
        sigpy_env(idx_range)  = y_seg(:,4);
        E_env(idx_range)      = y_seg(:,5);
        epsn2x_env(idx_range) = y_seg(:,6);
        epsn2y_env(idx_range) = y_seg(:,7);

        cur_idx = idx_range(end);
        current_s_slice = s_end;
    end

    if ode_failed
        continue;
    end

    all_slice_envelope(slice_idx).s = s_env(1:cur_idx);
    all_slice_envelope(slice_idx).sigx = sigx_env(1:cur_idx);
    all_slice_envelope(slice_idx).sigpx = sigpx_env(1:cur_idx);
    all_slice_envelope(slice_idx).sigy = sigy_env(1:cur_idx);
    all_slice_envelope(slice_idx).sigpy = sigpy_env(1:cur_idx);
    all_slice_envelope(slice_idx).E = E_env(1:cur_idx);
    all_slice_envelope(slice_idx).epsn2x = epsn2x_env(1:cur_idx);
    all_slice_envelope(slice_idx).epsn2y = epsn2y_env(1:cur_idx);
    all_slice_envelope(slice_idx).particle_count = slice_particles;
    all_slice_envelope(slice_idx).valid = true;

    [s_exit_uniq, u_idx] = unique(s_transport, 'stable');
    all_slice_centroid(slice_idx).s_exit  = s_exit_uniq;
    all_slice_centroid(slice_idx).x_exit  = xc(u_idx);
    all_slice_centroid(slice_idx).px_exit = pxc(u_idx);
    all_slice_centroid(slice_idx).y_exit  = yc(u_idx);
    all_slice_centroid(slice_idx).py_exit = pyc(u_idx);
    all_slice_centroid(slice_idx).particle_count = slice_particles;
    all_slice_centroid(slice_idx).valid = true;

    if isempty(s_envelope_common)
        s_envelope_common = s_env(1:cur_idx);
    end
end

% ================= 极速后处理与插值投影 =================
if isempty(s_envelope_common)
    s_envelope_common = s_transport;
end
num_dense = length(s_envelope_common);

% ================= 峰值流强切片几何发射度 (插值到 s_envelope_common) =================
peak_slice_emitx = zeros(num_dense, 1);
peak_slice_emity = zeros(num_dense, 1);

valid_mask_for_peak = false(num_slice, 1);
for ii = 1:num_slice
    valid_mask_for_peak(ii) = all_slice_envelope(ii).valid;
end
if any(valid_mask_for_peak)
    counts_peak = -inf(num_slice, 1);
    counts_peak(valid_mask_for_peak) = slice_particle_counts(valid_mask_for_peak);
    [~, peak_slice_idx] = max(counts_peak);

    s_pk     = all_slice_envelope(peak_slice_idx).s;
    gamma_pk = max(all_slice_envelope(peak_slice_idx).E, me_c2) / me_c2;
    emitx_pk = sqrt(max(all_slice_envelope(peak_slice_idx).epsn2x, 0)) ./ gamma_pk;
    emity_pk = sqrt(max(all_slice_envelope(peak_slice_idx).epsn2y, 0)) ./ gamma_pk;

    [s_pk_u, u_pk] = unique(s_pk, 'stable');
    peak_slice_emitx = interp1(s_pk_u, emitx_pk(u_pk), s_envelope_common, 'linear', 'extrap');
    peak_slice_emity = interp1(s_pk_u, emity_pk(u_pk), s_envelope_common, 'linear', 'extrap');
end

sum_w = zeros(num_dense, 1);
sum_x_c = zeros(num_dense, 1);
sum_y_c = zeros(num_dense, 1);
sum_px_c = zeros(num_dense, 1);
sum_py_c = zeros(num_dense, 1);
sum_gam = zeros(num_dense, 1);

for idx = 1:length(valid_slice_indices)
    slice_idx = valid_slice_indices(idx);
    if ~all_slice_centroid(slice_idx).valid
        continue;
    end
    w = all_slice_centroid(slice_idx).particle_count;

    x_i  = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).x_exit,  s_envelope_common, 'linear', 'extrap');
    y_i  = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).y_exit,  s_envelope_common, 'linear', 'extrap');
    px_i = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).px_exit, s_envelope_common, 'linear', 'extrap');
    py_i = interp1(all_slice_centroid(slice_idx).s_exit, all_slice_centroid(slice_idx).py_exit, s_envelope_common, 'linear', 'extrap');
    E_i  = interp1(all_slice_envelope(slice_idx).s,      all_slice_envelope(slice_idx).E,       s_envelope_common, 'linear', 'extrap');

    sum_w    = sum_w + w;
    sum_x_c  = sum_x_c + x_i*w;
    sum_y_c  = sum_y_c + y_i*w;
    sum_px_c = sum_px_c + px_i*w;
    sum_py_c = sum_py_c + py_i*w;
    sum_gam  = sum_gam + (E_i / me_c2)*w;
end

mean_x_c  = sum_x_c ./ sum_w;
mean_y_c  = sum_y_c ./ sum_w;
mean_px_c = sum_px_c ./ sum_w;
mean_py_c = sum_py_c ./ sum_w;
mean_gam  = sum_gam ./ sum_w;

sum_x2  = zeros(num_dense, 1);
sum_xxp = zeros(num_dense, 1);
sum_xp2 = zeros(num_dense, 1);
sum_y2  = zeros(num_dense, 1);
sum_yyp = zeros(num_dense, 1);
sum_yp2 = zeros(num_dense, 1);

for idx = 1:length(valid_slice_indices)
    slice_idx = valid_slice_indices(idx);
    if ~all_slice_envelope(slice_idx).valid
        continue;
    end
    w = all_slice_centroid(slice_idx).particle_count;

    s_c = all_slice_centroid(slice_idx).s_exit;
    s_e = all_slice_envelope(slice_idx).s;

    sx_i  = interp1(s_e, all_slice_envelope(slice_idx).sigx,   s_envelope_common, 'linear', 'extrap');
    spx_i = interp1(s_e, all_slice_envelope(slice_idx).sigpx,  s_envelope_common, 'linear', 'extrap');
    sy_i  = interp1(s_e, all_slice_envelope(slice_idx).sigy,   s_envelope_common, 'linear', 'extrap');
    spy_i = interp1(s_e, all_slice_envelope(slice_idx).sigpy,  s_envelope_common, 'linear', 'extrap');
    E_i   = interp1(s_e, all_slice_envelope(slice_idx).E,      s_envelope_common, 'linear', 'extrap');

    gam_v = E_i / me_c2;
    e2x_i = interp1(s_e, all_slice_envelope(slice_idx).epsn2x, s_envelope_common, 'linear', 'extrap');
    e2y_i = interp1(s_e, all_slice_envelope(slice_idx).epsn2y, s_envelope_common, 'linear', 'extrap');

    emx2 = max(e2x_i, 0) ./ (gam_v.^2);
    emy2 = max(e2y_i, 0) ./ (gam_v.^2);

    x_c  = interp1(s_c, all_slice_centroid(slice_idx).x_exit,  s_envelope_common, 'linear', 'extrap');
    px_c = interp1(s_c, all_slice_centroid(slice_idx).px_exit, s_envelope_common, 'linear', 'extrap');
    y_c  = interp1(s_c, all_slice_centroid(slice_idx).y_exit,  s_envelope_common, 'linear', 'extrap');
    py_c = interp1(s_c, all_slice_centroid(slice_idx).py_exit, s_envelope_common, 'linear', 'extrap');

    dx  =0;%x_c  - mean_x_c;
    dpx =0;%px_c - mean_px_c;
    dy  =0;%y_c  - mean_y_c;
    dpy =0;%py_c - mean_py_c;

    x2_s  = sx_i.^2 + dx.^2;
    y2_s  = sy_i.^2 + dy.^2;
    xxp_s = (sx_i.*spx_i) + dx.*dpx;
    yyp_s = (sy_i.*spy_i) + dy.*dpy;

    xp2_s = (emx2 + (sx_i.*spx_i).^2) ./ (sx_i.^2) + dpx.^2;
    yp2_s = (emy2 + (sy_i.*spy_i).^2) ./ (sy_i.^2) + dpy.^2;

    sum_x2  = sum_x2  + x2_s*w;
    sum_xxp = sum_xxp + xxp_s*w;
    sum_xp2 = sum_xp2 + xp2_s*w;
    sum_y2  = sum_y2  + y2_s*w;
    sum_yyp = sum_yyp + yyp_s*w;
    sum_yp2 = sum_yp2 + yp2_s*w;
end

Sigma_x  = sqrt(sum_x2 ./ sum_w);
Sigma_y  = sqrt(sum_y2 ./ sum_w);
Sigma_px = (sum_xxp ./ sum_w) ./ Sigma_x;
Sigma_py = (sum_yyp ./ sum_w) ./ Sigma_y;

Emit_x_proj = sqrt(max((sum_x2.*sum_xp2)./(sum_w.^2) - (sum_xxp./sum_w).^2, 0));
Emit_y_proj = sqrt(max((sum_y2.*sum_yp2)./(sum_w.^2) - (sum_yyp./sum_w).^2, 0));

results = struct( ...
    's_array', s_envelope_common, ...
    'sigx_array', Sigma_x, ...
    'sigy_array', Sigma_y, ...
    'sigpx_array', Sigma_px, ...
    'sigpy_array', Sigma_py, ...
    'x_centroid', mean_x_c, ...
    'y_centroid', mean_y_c, ...
    'px_centroid', mean_px_c, ...
    'py_centroid', mean_py_c, ...
    'gamma_array', mean_gam, ...
    'emitx_proj_array', Emit_x_proj, ...
    'emity_proj_array', Emit_y_proj, ...
    'k1_values', k1_values, ...
    'peak_slice_emitx', peak_slice_emitx, ...
    'peak_slice_emity', peak_slice_emity);
end

% ========================================================================
% ODE 核心引擎：空间电荷与包络积分
% ========================================================================
function dydt = envelopeODE_sliced(s, y, k1, A, phi, N_slice_real, z_val_in, rc, me_c2, ...
                                   s_start_elem, kz_drift, evolve_zval_flag)

    sigx = max(y(1), 1e-12);
    sigy = max(y(3), 1e-12);
    E    = max(y(5), me_c2);
    gamma = E / me_c2;
    epsn2x = max(y(6), 1e-30);
    epsn2y = max(y(7), 1e-30);

    if ~isfinite(sigx) || ~isfinite(sigy) || sigx <= 0 || sigy <= 0
        dydt = zeros(7,1);
        return;
    end

    emitx2_geo = epsn2x / (gamma^2);
    emittance_x = emitx2_geo / (sigx^3);

    emity2_geo = epsn2y / (gamma^2);
    emittance_y = emity2_geo / (sigy^3);

    % --- 只让 z_val 进入 Ksc ---
    if evolve_zval_flag
        z_val_local = z_val_in * localZScale(s - s_start_elem, kz_drift);
    else
        z_val_local = z_val_in;
    end

    sc_x = 0;
    sc_y = 0;
    if N_slice_real > 0 && z_val_local > 0
        lambda_z = N_slice_real / z_val_local;
        K_perv   = (2 * rc * lambda_z) / (gamma^3);
        sc_x     = K_perv / (2* 1*(sigx + sigy));
        sc_y     = K_perv / (2 *1.* (sigx + sigy));
    end

    if A ~= 0
        dE_ds = A * cos(phi);
        dgamma_over_gamma = dE_ds / E;
        second_order_coeff = - (A / E)^2 /8;
    else
        dE_ds = 0;
        dgamma_over_gamma = 0;
        second_order_coeff = 0;
    end

    dsigx_ds  = y(2);
    dsigy_ds  = y(4);
    dsigpx_ds = 1*emittance_x + sc_x - k1*sigx - dgamma_over_gamma*y(2) + second_order_coeff*sigx;
    dsigpy_ds = 1*emittance_y + sc_y + k1*sigy - dgamma_over_gamma*y(4) + second_order_coeff*sigy;

    % --- 切片归一化发射度平方增长 (Wangler-Crandall-Mills) ---
    % 在 ODE 内部以连续速率累加，使 emittance_x/y 在同一元件内自洽更新
    depsn2x_ds = 0;
    depsn2y_ds = 0;
    if N_slice_real > 0 && z_val_local > 0
        ds_kick_local = 0.02;   % 与 odeset MaxStep 同步
        A0 = (2 * rc * (N_slice_real / z_val_local)) / (gamma^3);

        kxy = min(max(sigy/sigx, 1e-6), 1e6);
        kyx = min(max(sigx/sigy, 1e-6), 1e6);

        ux = 0.5 * min(max(kxy^(-1/5.5), 1e-12), 1e12);
        uy = 0.5 * min(max(kyx^(-1/5.5), 1e-12), 1e12);
        Fx = (2*sqrt(1+2*ux) - sqrt(1+4*ux) - 1) - (1 - 1/sqrt(1+2*ux))^2;
        Fy = (2*sqrt(1+2*uy) - sqrt(1+4*uy) - 1) - (1 - 1/sqrt(1+2*uy))^2;
        Fx = max(Fx, 0); Fy = max(Fy, 0);

        exp_arg_x = max(-1/(8*kxy^(1/5.5)), -50);
        g_x = (1/4.25) / (1 + kxy) / (1 - max(exp(exp_arg_x), 1e-15));
        enxsc = (g_x * A0)^2;

        exp_arg_y = max(-1/(8*kyx^(1/5.5)), -50);
        g_y = (1/4.25) / (1 + kyx) / (1 - max(exp(exp_arg_y), 1e-15));
        enysc = (g_y * A0)^2;

        depsn2x_ds = 1*gamma^2 * enxsc * Fx * ds_kick_local;
        depsn2y_ds = 1*gamma^2 * enysc * Fy * ds_kick_local;
    end

    dydt = [dsigx_ds; dsigpx_ds; dsigy_ds; dsigpy_ds; dE_ds; depsn2x_ds; depsn2y_ds];
end

% ========================================================================
% (legacy) post-processing growth — superseded by in-ODE 增长，保留供调试
% ========================================================================
function [eps2x_new, eps2y_new] = applySliceEmitGrowth(z_seg, y_seg, N_slice_real, z_val_in, rc, me_c2) %#ok<DEFNU>
    eps2x_new = y_seg(:, 6);
    eps2y_new = y_seg(:, 7);
    n = length(z_seg);
    if n < 2 || N_slice_real <= 0 || z_val_in <= 0
        return;
    end

    sigx_arr  = max(y_seg(:,1), 1e-12);
    sigy_arr  = max(y_seg(:,3), 1e-12);
    gamma_arr = max(y_seg(:,5), me_c2) / me_c2;
    A0_arr    = (2 * rc * (N_slice_real / z_val_in)) ./ (gamma_arr.^3);

    for i = 1:(n-1)
        ds_local = z_seg(i+1) - z_seg(i);
        if ds_local <= 0, continue; end

        sigx  = sigx_arr(i);  sigy  = sigy_arr(i);
        gamma = gamma_arr(i); A0    = A0_arr(i);

        kxy = min(max(sigy/sigx, 1e-6), 1e6);
        kyx = min(max(sigx/sigy, 1e-6), 1e6);

        ux = 0.5 * min(max(kxy^(-1/5.5), 1e-12), 1e12);
        uy = 0.5 * min(max(kyx^(-1/5.5), 1e-12), 1e12);
        Fx = (2*sqrt(1+2*ux) - sqrt(1+4*ux) - 1) - (1 - 1/sqrt(1+2*ux))^2;
        Fy = (2*sqrt(1+2*uy) - sqrt(1+4*uy) - 1) - (1 - 1/sqrt(1+2*uy))^2;
        Fx = max(Fx, 0); Fy = max(Fy, 0);

        exp_arg_x = max(-1/(8*kxy^(1/5.5)), -50);
        g_x = (1/4.25) / (1 + kxy) / (1 - max(exp(exp_arg_x), 1e-15));
        enxsc = (g_x * A0)^2;

        exp_arg_y = max(-1/(8*kyx^(1/5.5)), -50);
        g_y = (1/4.25) / (1 + kyx) / (1 - max(exp(exp_arg_y), 1e-15));
        enysc = (g_y * A0)^2;

        eps2x_new(i+1) = eps2x_new(i) + gamma^2 * enxsc * Fx * ds_local/6;
        eps2y_new(i+1) = eps2y_new(i) + gamma^2 * enysc * Fy * ds_local/6;
    end
end

% ========================================================================
% 低能 drift 的速度压缩等效 R56（首版近似）
% 注意：符号和你的 z 定义有关，必要时把 sign_flag 改号 benchmark
% ========================================================================
function R56v = driftR56_velocity(L, gamma)
    beta = sqrt(max(1 - 1/gamma^2, 1e-14));
    R56v = L / (beta^2 * gamma^2);
end

% ========================================================================
% z_val 的局部线性缩放：z_val(s) = z_val_in * |1 + kz * s|
% ========================================================================
function scale = localZScale(ds_local, kz)
    scale = abs(1 + kz * ds_local);
    scale = max(scale, 1e-6);
end

function val = spaceChargeXFast(sigx, sigy, A0, sigx_inv)
    k = min(max(sigy*sigx_inv, 1e-6), 1e6);
    k_power = min(max(k^(-1/5.5), 1e-6), 1e6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    val = min(max(f * (1/4.25)/(1+k)/(1-max(exp(exp_arg), 1e-15)) * A0 * sigx_inv, -1e6), 1e6);
end

function val = spaceChargeYFast(sigx, sigy, A0, sigy_inv)
    k = min(max(sigx*sigy_inv, 1e-6), 1e6);
    k_power = min(max(k^(-1/5.5), 1e-6), 1e6);
    f = 1 - 1/sqrt(1 + k_power);
    exp_arg = max(-1/(8*k^(1/5.5)), -50);
    val = min(max(f * (1/4.25)/(1+k)/(1-max(exp(exp_arg), 1e-15)) * A0 * sigy_inv, -1e6), 1e6);
end