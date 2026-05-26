% 清除工作区和关闭所有图表
clear all;
close all;

% 基本参数
N = 1e10;                     % 粒子数
e = 1.602e-19;                % 电子电荷 (C)
epsilon_0 = 8.85e-12;         % 真空介电常数 (F/m)
sigma_z = 1e-3;               % z方向标准差 (m)

% 创建z坐标范围
z_range = linspace(-5*sigma_z, 5*sigma_z, 200);

% 不同的sigmax/sigmay比值
ratio_values = [1/30, 5, 15, 30];
n_ratios = length(ratio_values);

% 预分配场强数组
Ex_values = zeros(length(z_range), n_ratios);

% 创建从红[.6,0,0]到蓝[0,0,1]的渐变颜色
colors = zeros(n_ratios, 3);
for i = 1:n_ratios
    t = (i-1)/(n_ratios-1);  % 从0到1的归一化参数
    colors(i,:) = [.6*(1-t), 0, t];  % 从[.6,0,0]渐变到[0,0,1]
end

% 积分函数 - 电势
potential_integrand = @(q, x, y, z, sigma_x, sigma_y, sigma_z) ...
    exp(-x^2/(2*sigma_x^2+q) - y^2/(2*sigma_y^2+q) - z^2/(2*sigma_z^2+q)) ./ ...
    sqrt((2*sigma_x^2+q).*(2*sigma_y^2+q).*(2*sigma_z^2+q));

% 对各种比值和z值计算场强
for ratio_idx = 1:n_ratios
    ratio = ratio_values(ratio_idx);
    sigma_x = 1e-3;              % x方向标准差 (m)
    sigma_y = sigma_x / ratio;   % y方向标准差 (m)
    
    fprintf('计算 sigma_x/sigma_y = %g 的场强\n', ratio);
    
    % 计算每个z位置的场强
    for i = 1:length(z_range)
        z = z_range(i);
        x = sigma_x;  % 在r=sigma_x处计算Ex
        y = 0;        % y=0平面
        
        % 计算电势的数值积分
        q_max = 20*max([sigma_x^2, sigma_y^2, sigma_z^2]);  % 积分上限
        q_points = 1000;
        q_values = linspace(0, q_max, q_points);
        dq = q_values(2) - q_values(1);
        
        % 数值积分计算电势
        potential_sum = 0;
        for q_idx = 1:length(q_values)
            q = q_values(q_idx);
            integrand = potential_integrand(q, x, y, z, sigma_x, sigma_y, sigma_z);
            potential_sum = potential_sum + integrand * dq;
        end
        
        % 计算电势
        potential = N * e / (4 * pi * epsilon_0 * sqrt(pi)) * potential_sum;
        
        % 计算Ex，使用中心差分近似导数
        dx = sigma_x * 0.01;
        x_plus = x + dx;
        x_minus = x - dx;
        
        % 计算x+dx处的电势
        potential_plus_sum = 0;
        for q_idx = 1:length(q_values)
            q = q_values(q_idx);
            integrand = potential_integrand(q, x_plus, y, z, sigma_x, sigma_y, sigma_z);
            potential_plus_sum = potential_plus_sum + integrand * dq;
        end
        potential_plus = N * e / (4 * pi * epsilon_0 * sqrt(pi)) * potential_plus_sum;
        
        % 计算x-dx处的电势
        potential_minus_sum = 0;
        for q_idx = 1:length(q_values)
            q = q_values(q_idx);
            integrand = potential_integrand(q, x_minus, y, z, sigma_x, sigma_y, sigma_z);
            potential_minus_sum = potential_minus_sum + integrand * dq;
        end
        potential_minus = N * e / (4 * pi * epsilon_0 * sqrt(pi)) * potential_minus_sum;
        
        % 使用中心差分计算Ex
        Ex = -(potential_plus - potential_minus) / (2 * dx);
        
        % 存储结果
        Ex_values(i, ratio_idx) = Ex;
        
        % 显示进度
        if mod(i, 20) == 0
            fprintf('已计算 %d/%d 个点\n', i, length(z_range));
        end
    end
end

% 创建高斯分布 - 归一化到每个场分布的最大值
gaussian = @(z, sigma_z) exp(-(z.^2)/(2*sigma_z^2));
gaussian_values = gaussian(z_range, sigma_z);

% 绘制结果 - 使用每个曲线自己的最大值归一化
figure('Position', [100, 100, 900, 600]);
hold on;

% 先画归一化后的高斯分布（虚线）
plot(z_range/sigma_z, gaussian_values, '--k', 'LineWidth', 1.5);

% 逐个画归一化的场强分布
for ratio_idx = 1:n_ratios
    % 找到当前比值的最大场强
    max_Ex = max(abs(Ex_values(:, ratio_idx)));
    normalized_Ex = Ex_values(:, ratio_idx) / max_Ex;
    
    % 绘制归一化场强
    plot(z_range/sigma_z, normalized_Ex, 'LineWidth', 2, 'Color', colors(ratio_idx,:));
end

hold off;
grid on;
xlabel('z/\sigma_z', 'FontName', 'Cambria Math', 'FontSize', 16);
ylabel('归一化Ex场强 (在r=\sigma_x处)', 'FontName', 'Cambria Math', 'FontSize', 16);
title('不同横向纵横比下的纵向电场分布', 'FontName', 'Cambria Math', 'FontSize', 18);

% 添加图例
legend_strs = cell(n_ratios + 1, 1);
legend_strs{1} = '高斯分布';
for ratio_idx = 1:n_ratios
    ratio = ratio_values(ratio_idx);
    legend_strs{ratio_idx + 1} = ['\sigma_x/\sigma_y = ', num2str(ratio)];
end
legend(legend_strs, 'Location', 'best', 'FontName', 'Cambria Math', 'FontSize', 12);
set(gca, 'FontName', 'Cambria Math', 'FontSize', 14);

% 保存结果
save('Ex_field_results.mat', 'z_range', 'sigma_z', 'Ex_values', 'ratio_values');