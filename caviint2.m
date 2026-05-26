clear; clc; close all;

% ==========================================
% 1. 物理参数定义 (Physical Parameters)
% ==========================================
L = 1.0;                % 腔体长度 (m)
phi_deg = 0;           % RF相位 (度)
Delta_gamma = 20;       % 能量增益
gamma_start = 1.5;      % 初始扫描能量
gamma_end = 20.0;       % 终止扫描能量

% 计算常数
phi = deg2rad(phi_deg);
Omega = 1 / (sqrt(8) * cos(phi));
K = sqrt(2) * cos(phi);

% ==========================================
% 2. 开始扫描计算
% ==========================================
gamma_in_list = linspace(gamma_start, gamma_end, 50);

% 预分配内存
I1_exact = zeros(size(gamma_in_list));
I1_0ord  = zeros(size(gamma_in_list));
I2_exact = zeros(size(gamma_in_list));
I2_0ord  = zeros(size(gamma_in_list));

fprintf('正在进行数值积分与 0阶解析公式 对比...\n');

for k = 1:length(gamma_in_list)
    gi = gamma_in_list(k);
    gf = gi + Delta_gamma;
    gamma_prime = (gf - gi) / L;
    
    % 入口处的总相位推进
    alpha_i = Omega * log(gf / gi);
    
    % ---------------------------------------------------
    % 1. 计算 I1 (对应 R12)
    % ---------------------------------------------------
    % [数值精确积分]
    % 注意被积函数包含 1/sqrt(g^2 - 1)
    C1 = 1 / (Omega * gamma_prime^2 * sqrt(gf));
    func_I1 = @(g) C1 .* (sqrt(g) .* sin(Omega .* log(gf ./ g))) ./ sqrt(g.^2 - 1);
    I1_exact(k) = integral(func_I1, gi, gf);
    
    % [0阶解析代数公式]
    term1_0ord = 1 / (gamma_prime^2 * (0.25 + Omega^2));
    term2_0ord = 1 - sqrt(gi/gf) * (cos(alpha_i) + (1/(2*Omega)) * sin(alpha_i));
    I1_0ord(k) = term1_0ord * term2_0ord;
    
    % ---------------------------------------------------
    % 2. 计算 I2 (对应 R22)
    % ---------------------------------------------------
    % [数值精确积分] -> 修正：C2 分母只能有一次方 gamma_prime
    C2 = 1 / (gamma_prime * gf);
    func_I2 = @(g) C2 .* (g .* (cos(Omega .* log(gf ./ g)) + K .* sin(Omega .* log(gf ./ g)))) ./ sqrt(g.^2 - 1);
    I2_exact(k) = integral(func_I2, gi, gf);
    
    % [0阶解析代数公式] -> 修正：term1_I2_0ord 分母只能有一次方 gamma_prime
    term1_I2_0ord = 1 / (gamma_prime * (1 + Omega^2));
    term2_I2_0ord = 1.5 - (gi/gf) * (1.5 * cos(alpha_i) + (K - Omega) * sin(alpha_i));
    I2_0ord(k) = term1_I2_0ord * term2_I2_0ord;
end

% 计算相对误差 (%)
err_I1 = abs(I1_0ord - I1_exact) ./ abs(I1_exact) * 100;
err_I2 = abs(I2_0ord - I2_exact) ./ abs(I2_exact) * 100;

% ==========================================
% 3. 绘图 (Visualization)
% ==========================================
figure('Color', 'w', 'Position', [100, 100, 1000, 600]);

% --- 子图 1: I1 积分值对比 ---
subplot(2, 2, 1);
plot(gamma_in_list, I1_exact, 'k-', 'LineWidth', 2); hold on;
plot(gamma_in_list, I1_0ord, 'r--', 'LineWidth', 2);
title('I_1 (R_{12}) 积分对比');
xlabel('\gamma_{in}'); ylabel('积分值');
legend('精确数值积分', '0阶代数公式', 'Location', 'best');
grid on;

% --- 子图 2: I1 相对误差 ---
subplot(2, 2, 2);
semilogy(gamma_in_list, err_I1, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
yline(1, 'k:', '1% 误差线');
title('I_1 (R_{12}) 相对误差 (%)');
xlabel('\gamma_{in}'); ylabel('误差 (%)');
grid on;

% --- 子图 3: I2 积分值对比 ---
subplot(2, 2, 3);
plot(gamma_in_list, I2_exact, 'k-', 'LineWidth', 2); hold on;
plot(gamma_in_list, I2_0ord, 'b--', 'LineWidth', 2);
title('I_2 (R_{22}) 积分对比');
xlabel('\gamma_{in}'); ylabel('积分值');
legend('精确数值积分', '0阶代数公式', 'Location', 'best');
grid on;

% --- 子图 4: I2 相对误差 ---
subplot(2, 2, 4);
semilogy(gamma_in_list, err_I2, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 4);
yline(1, 'k:', '1% 误差线');
title('I_2 (R_{22}) 相对误差 (%)');
xlabel('\gamma_{in}'); ylabel('误差 (%)');
grid on;

% 打印关键点的误差信息
fprintf('\n在极低能端 (gamma_in = 1.5) 的表现：\n');
fprintf('  I1 (R12) 0阶公式相对误差: %.2f%%\n', err_I1(1));
fprintf('  I2 (R22) 0阶公式相对误差: %.2f%%\n', err_I2(1));

fprintf('\n在常规能端 (gamma_in = 10.0) 的表现：\n');
idx_10 = find(gamma_in_list >= 10, 1);
fprintf('  I1 (R12) 0阶公式相对误差: %.2f%%\n', err_I1(idx_10));
fprintf('  I2 (R22) 0阶公式相对误差: %.2f%%\n', err_I2(idx_10));