function slice_data = read_z_distribution(infile, num_slices, threshold_n, n_sigma)
    % ---------------------------------------------------------------------
    % read_z_distribution  —  读取 beam 文件，3sigma 截断，切片分析
    % ---------------------------------------------------------------------
    % 默认值设置
    if nargin < 2 || isempty(num_slices),  num_slices  = 20;  end
    if nargin < 3 || isempty(threshold_n), threshold_n = 0;   end
    if nargin < 4 || isempty(n_sigma),     n_sigma     = 3;   end
    
    % 防御性编程：防止通过命令行调用时参数被误认为字符串
    if ischar(num_slices)  || isstring(num_slices),  num_slices  = str2double(num_slices);  end
    if ischar(threshold_n) || isstring(threshold_n), threshold_n = str2double(threshold_n); end
    if ischar(n_sigma)     || isstring(n_sigma),     n_sigma     = str2double(n_sigma);     end
    
    % =====================================================================
    % 1) 读取文件
    % =====================================================================
    [np_orig] = textread(infile,'%f%*s%*s%*s',1,'headerlines',3);
    [Q_bunch] = textread(infile,'%f%*s%*s%*s',1,'headerlines',6);
    [x, xp, y, yp, z, p] = textread(infile, ...
        '%f%f%f%f%f%f%*f%*f%*f%*f%*f%*f%*f', np_orig, 'headerlines', 9);
    
    % 确保读取出来的一定是列向量
    x = x(:); xp = xp(:); y = y(:); yp = yp(:); z = z(:); p = p(:);
    delta = p;
    
    fprintf('Total particles read = %d\n', numel(x));
    fprintf('Bunch charge from file = %.6e C\n', Q_bunch);
    
    % =====================================================================
    % 2) 横向 n_sigma 截断 (基于全束团 rms)
    % =====================================================================
    x_mean = mean(x);  y_mean = mean(y);
    sig_x  = std(x,1);  sig_y  = std(y,1);
    
    keep_xy = (abs(x - x_mean) <= n_sigma * sig_x) & ...
              (abs(y - y_mean) <= n_sigma * sig_y);
              
    np_before_cut = numel(x);
    x     = x(keep_xy);
    xp    = xp(keep_xy);
    y     = y(keep_xy);
    yp    = yp(keep_xy);
    z     = z(keep_xy);
    delta = delta(keep_xy);
    
    np_after_cut = numel(x);
    
    % 截断后重新计算 rms
    sig_x_cut = std(x,1);
    sig_y_cut = std(y,1);
    
    fprintf('\n--- Transverse %d-sigma cut ---\n', n_sigma);
    fprintf('  Before: %d particles,  sigma_x = %.6e,  sigma_y = %.6e\n', ...
        np_before_cut, sig_x, sig_y);
    fprintf('  After : %d particles,  sigma_x = %.6e,  sigma_y = %.6e\n', ...
        np_after_cut, sig_x_cut, sig_y_cut);
    fprintf('  Removed: %d particles (%.2f%%)\n', ...
        np_before_cut - np_after_cut, 100*(1 - np_after_cut/np_before_cut));
        
    % 截断后的有效电荷
    Q_cut = Q_bunch * np_after_cut / np_orig;
    fprintf('  Effective charge after cut = %.6e C\n', Q_cut);
    
    % =====================================================================
    % 3) 参考能量 E0 (中心切片 delta=0 附近的平均能量)
    % =====================================================================
    me_c2 = 0.511e6;  % eV
    [~, idx_ref] = min(abs(delta));
    E0_eV = me_c2 * (1 + delta(idx_ref));
    
    z_center = mean(z);
    near_center = abs(z - z_center) < std(z,1) * 0.2;
    if sum(near_center) > 10
        E0_eV = me_c2 * (1 + mean(delta(near_center)));
    end
    fprintf('  Reference energy E0 = %.6e eV\n', E0_eV);
    
    % =====================================================================
    % 4) 纵向切分
    % =====================================================================
    z_min = min(z);
    z_max = max(z);
    z_length_orig = z_max - z_min;
    dz = z_length_orig / num_slices;
    z_edges = linspace(z_min, z_max, num_slices + 1);
    
    [N_counts, ~, bin_idx] = histcounts(z, z_edges);
    
    % =====================================================================
    % 5) 流强截断
    % =====================================================================
    N_max = max(N_counts);
    thresh_N = threshold_n * N_max;
    
    valid_bins = find(N_counts >= thresh_N);
    num_valid = length(valid_bins);
    
    core_mask = ismember(bin_idx, valid_bins);
    np_core = sum(core_mask);
    valid_z_length = num_valid * dz;
    
    if threshold_n > 0
        fprintf('\n--- Longitudinal current threshold (n = %.2f * peak) ---\n', threshold_n);
        fprintf('  Total slices: %d,  Kept: %d,  Removed: %d\n', ...
            num_slices, num_valid, num_slices - num_valid);
        fprintf('  Valid Z length: %.6e m\n', valid_z_length);
    end
    
    % =====================================================================
    % 6) 组装结果矩阵
    % =====================================================================
    slice_data = zeros(num_valid + 2, 17);
    
    % --- 第一行：全局统计 ---
    slice_data(1, 1)  = np_after_cut;
    slice_data(1, 10) = mean(x(core_mask));
    slice_data(1, 11) = mean(xp(core_mask));
    slice_data(1, 12) = mean(y(core_mask));
    slice_data(1, 13) = mean(yp(core_mask));
    slice_data(1, 14) = mean(z(core_mask));
    slice_data(1, 15) = mean(delta(core_mask));
    slice_data(1, 16) = Q_cut;
    slice_data(1, 17) = E0_eV;
    
    % --- 最后一行 ---
    slice_data(end, 1) = valid_z_length;
    
    % --- 切片行 ---
    for k = 1:num_valid
        b = valid_bins(k);
        slice_idx = (bin_idx == b);
        slice_particles = sum(slice_idx);
        slice_data(k+1, 1) = slice_particles;
        
        if slice_particles > 0
            x_slice     = x(slice_idx);
            xp_slice    = xp(slice_idx);
            y_slice     = y(slice_idx);
            yp_slice    = yp(slice_idx);
            z_slice     = z(slice_idx);
            delta_slice = delta(slice_idx);
            
            x_mean_s  = mean(x_slice);
            xp_mean_s = mean(xp_slice);
            y_mean_s  = mean(y_slice);
            yp_mean_s = mean(yp_slice);
            z_mean_s  = mean(z_slice);
            delta_mean_s = mean(delta_slice);
            
            slice_data(k+1, 10) = x_mean_s;
            slice_data(k+1, 11) = xp_mean_s;
            slice_data(k+1, 12) = y_mean_s;
            slice_data(k+1, 13) = yp_mean_s;
            slice_data(k+1, 14) = z_mean_s;
            slice_data(k+1, 15) = delta_mean_s;
            
            slice_data(k+1, 16) = mean(x_slice .* delta_slice);
            slice_data(k+1, 17) = mean(xp_slice .* delta_slice);
            
            x_rel  = x_slice - x_mean_s;
            xp_rel = xp_slice - xp_mean_s;
            y_rel  = y_slice - y_mean_s;
            yp_rel = yp_slice - yp_mean_s;
            
            sigma11_x = mean(x_rel.^2);
            sigma12_x = mean(x_rel .* xp_rel);
            sigma22_x = mean(xp_rel.^2);
            sigma11_y = mean(y_rel.^2);
            sigma12_y = mean(y_rel .* yp_rel);
            sigma22_y = mean(yp_rel.^2);
            
            emit_x_geo = sqrt(max(0, sigma11_x * sigma22_x - sigma12_x^2));
            emit_y_geo = sqrt(max(0, sigma11_y * sigma22_y - sigma12_y^2));
            
            if emit_x_geo < 1e-15, emit_x_geo = 1e-15; end
            if emit_y_geo < 1e-15, emit_y_geo = 1e-15; end
            
            slice_data(k+1, 2) = sigma11_x / emit_x_geo;
            slice_data(k+1, 3) = -sigma12_x / emit_x_geo;
            slice_data(k+1, 4) = sigma22_x / emit_x_geo;
            slice_data(k+1, 5) = emit_x_geo;
            slice_data(k+1, 6) = sigma11_y / emit_y_geo;
            slice_data(k+1, 7) = -sigma12_y / emit_y_geo;
            slice_data(k+1, 8) = sigma22_y / emit_y_geo;
            slice_data(k+1, 9) = emit_y_geo;
        else
            slice_data(k+1, 2:9) = [1.0, 0.0, 1.0, 1e-9, 1.0, 0.0, 1.0, 1e-9];
        end
    end
    
    % =====================================================================
    % 7) 保存文件
    % =====================================================================
    [filepath, name, ~] = fileparts(infile);
    output_filename = fullfile(filepath, [name, '_slice_params.txt']);
    
    fid = fopen(output_filename, 'w');
    if fid == -1
        error('Cannot create output file: %s', output_filename);
    end
    
    fprintf(fid, '%% Slice Optics Analysis Report\n');
    fprintf(fid, '%% Source: %s\n', infile);
    fprintf(fid, '%% Date: %s\n', datestr(now));
    fprintf(fid, '%% Transverse cut: %d sigma (removed %d / %d particles)\n', ...
        n_sigma, np_before_cut - np_after_cut, np_before_cut);
    fprintf(fid, '%% Slices: %d total, %d kept (threshold = %.2f)\n', ...
        num_slices, num_valid, threshold_n);
    fprintf(fid, '%% Slice thickness dz = %.6e m\n', dz);
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Row 1: global [col1=np, col16=Q_bunch(C), col17=E0(eV)]\n');
    fprintf(fid, '%% Row 2~%d: slice data\n', num_valid+1);
    fprintf(fid, '%% Last row: col1 = valid_z_length\n');
    fprintf(fid, '%% \n');
    fprintf(fid, '%% Col01:N  Col02:bx  Col03:ax  Col04:gx  Col05:ex  Col06:by  Col07:ay  Col08:gy  Col09:ey\n');
    fprintf(fid, '%% Col10:xc Col11:xpc Col12:yc Col13:ypc Col14:zc Col15:dc Col16:x_p/Q Col17:xp_p/E0\n');
    for i = 1:size(slice_data, 1)
        fprintf(fid, '%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n', slice_data(i, :));
    end
    fclose(fid);
    
    fprintf('\nSaved to: %s\n', output_filename);
end