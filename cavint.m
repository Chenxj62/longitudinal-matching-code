clear; clc; close all;

% ==========================================
% 1. Physics Parameters
% ==========================================
L = 1.0;              % Cavity length (m)
phi_deg = 20;         % RF phase (degrees)
Delta_gamma = 10;     % Energy gain
gamma_start = 2.0;    % Minimum input gamma (must be >= 2 for expansion)
gamma_end = 20.0;     % Maximum input gamma

% Derived parameters
phi = deg2rad(phi_deg);
% Focusing parameter Omega (Rosenzweig-Serafini model)
Omega = 1 / (sqrt(8) * cos(phi)); 

% ==========================================
% 2. Define Analytical Helper Functions
% ==========================================
% The integral form is: Int [ gamma^n * sin(Omega * ln(gf/g)) ] dg
% Solution: F(g) = (gamma^(n+1) / ((n+1)^2 + Omega^2)) * [ (n+1)*sin(alpha) + Omega*cos(alpha) ]
% We need to evaluate [F(gamma_f) - F(gamma_i)]

% Function to calculate the primitives
get_primitive = @(g, gf, n, Omg) (g.^(n+1) ./ ((n+1)^2 + Omg^2)) .* ...
    ( (n+1).*sin(Omg.*log(gf./g)) + Omg.*cos(Omg.*log(gf./g)) );

% ==========================================
% 3. Simulation Loop
% ==========================================
gamma_in_list = linspace(gamma_start, gamma_end, 50);

% Pre-allocate arrays
I_exact     = zeros(size(gamma_in_list));
I_standard  = zeros(size(gamma_in_list));
I_corrected = zeros(size(gamma_in_list));
err_standard = zeros(size(gamma_in_list));
err_corrected = zeros(size(gamma_in_list));

fprintf('Running simulation for gamma_in = [%.1f, %.1f]...\n', gamma_start, gamma_end);

for k = 1:length(gamma_in_list)
    g_i = gamma_in_list(k);
    g_f = g_i + Delta_gamma;
    gamma_prime = (g_f - g_i) / L;

    % Common prefactor for R12: 1 / (Omega * gamma' * sqrt(gamma_f))
    % Note: R12 itself has sqrt(gamma/gamma_f). We factor out constants.
    K_factor = 1 / (Omega * gamma_prime * sqrt(g_f));

    % --- Method A: Exact Numerical Integration ---
    % Integrand: R12 / sqrt(g^2 - 1)
    % R12 part: sqrt(g) * sin(Omega * ln(gf/g))
    func_exact = @(g) (1./sqrt(g.^2 - 1)) .* sqrt(g) .* sin(Omega .* log(g_f./g));
    I_exact(k) = K_factor * integral(func_exact, g_i, g_f);

    % --- Method B: Standard Analytical (Ultra-relativistic) ---
    % Approx: 1/sqrt(g^2-1) ~ 1/gamma
    % Effective Integrand power: g^(-1) * g^(0.5) = g^(-0.5) -> n = -0.5
    val_std_f = get_primitive(g_f, g_f, -0.5, Omega);
    val_std_i = get_primitive(g_i, g_f, -0.5, Omega);
    I_standard(k) = K_factor * (val_std_f - val_std_i);

    % --- Method C: Corrected Analytical (First Order Expansion) ---
    % Expansion: 1/sqrt(g^2-1) ~ 1/gamma + 1/(2*gamma^3)
    % We already have the first term (Method B). 
    % We add the second term: 0.5 * g^(-3) * g^(0.5) = 0.5 * g^(-2.5) -> n = -2.5

    val_corr_f = get_primitive(g_f, g_f, -2.5, Omega);
    val_corr_i = get_primitive(g_i, g_f, -2.5, Omega);

    term_correction = 0.5 * (val_corr_f - val_corr_i);
    I_corrected(k) = I_standard(k) + K_factor * term_correction;

    % --- Calculate Errors ---
    err_standard(k)  = abs(I_standard(k) - I_exact(k))  / abs(I_exact(k)) * 100;
    err_corrected(k) = abs(I_corrected(k) - I_exact(k)) / abs(I_exact(k)) * 100;
end

% ==========================================
% 4. Visualization
% ==========================================
figure('Color', 'w', 'Position', [100, 100, 1000, 400]);

% Plot 1: Integral Values
subplot(1, 2, 1);
plot(gamma_in_list, I_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Numerical'); hold on;
plot(gamma_in_list, I_standard, 'b--', 'LineWidth', 1.5, 'DisplayName', '1st ord');
plot(gamma_in_list, I_corrected, 'r-.', 'LineWidth', 1.5, 'DisplayName', '2nd ord');
xlabel('Input Energy \gamma_{in}');
ylabel('Integral Value (m/rad)');
title('\int R_{12}^{s->s_f} ds');
legend('Location', 'best');
grid on;

% Plot 2: Relative Error (Log Scale)
subplot(1, 2, 2);
semilogy(gamma_in_list, err_standard, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
semilogy(gamma_in_list, err_corrected, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);


xlabel('Input Energy \gamma_{in}');
ylabel('Relative Error (%)');
%title('Accuracy Improvement: Standard vs Corrected');
legend('1st ord Error', '2nd ord Error', 'Location', 'NorthEast');
grid on;



% clear; clc; close all;
% 
% % ==========================================
% % 1. 物理参数
% % ==========================================
% L = 1.0;
% phi_deg = 10;
% Delta_gamma = 10;
% gamma_start = 2.0;   % 入口能量
% gamma_end = 20.0;
% 
% phi = deg2rad(phi_deg);
% Omega = 1 / (sqrt(8) * cos(phi));
% % R22 中的正弦项系数 K (Rosenzweig-Serafini)
% K_sin = sqrt(2) * cos(phi);
% 
% % ==========================================
% % 2. 定义积分原函数 (Primitive Functions)
% % ==========================================
% % 通用积分公式: Int g^n * Trig(Omega * ln(gf/g)) dg
% % 这些是解析解的核心
% 
% % 对应 sin(alpha) 的积分
% get_prim_sin = @(g, gf, n, Omg) (g.^(n+1) ./ ((n+1)^2 + Omg^2)) .* ...
%     ( (n+1).*sin(Omg.*log(gf./g)) + Omg.*cos(Omg.*log(gf./g)) );
% 
% % 对应 cos(alpha) 的积分
% get_prim_cos = @(g, gf, n, Omg) (g.^(n+1) ./ ((n+1)^2 + Omg^2)) .* ...
%     ( (n+1).*cos(Omg.*log(gf./g)) - Omg.*sin(Omg.*log(gf./g)) );
% 
% % 组合: 计算 R22 核心部分 [cos + K*sin] 的积分
% get_R22_int = @(g, gf, n, Omg, K) get_prim_cos(g, gf, n, Omg) + K * get_prim_sin(g, gf, n, Omg);
% 
% % ==========================================
% % 3. 计算循环
% % ==========================================
% gamma_in_list = linspace(gamma_start, gamma_end, 50);
% I_exact = zeros(size(gamma_in_list));
% I_standard = zeros(size(gamma_in_list));
% I_corrected = zeros(size(gamma_in_list));
% 
% % 误差记录
% err_std = zeros(size(gamma_in_list));
% err_cor = zeros(size(gamma_in_list));
% 
% for k = 1:length(gamma_in_list)
%     g_i = gamma_in_list(k);
%     g_f = g_i + Delta_gamma;
%     gamma_prime = (g_f - g_i) / L;
% 
%     % R22 积分前的常数因子
%     % 原始积分: Int (1/gamma') * (R22 / sqrt(g^2-1)) dg
%     % R22 = (g/gf) * [Trig]
%     % 提取常数 C = 1 / (gamma' * gf)
%     C_factor = 1 / (gamma_prime * g_f);
% 
%     % --- A. 精确数值积分 (The Truth) ---
%     % 被积: g * [Trig] / sqrt(g^2-1)
%     alpha_func = @(g) Omega .* log(g_f ./ g);
%     trig_part  = @(g) cos(alpha_func(g)) + K_sin * sin(alpha_func(g));
%     func_exact = @(g) (g ./ sqrt(g.^2 - 1)) .* trig_part(g);
% 
%     I_exact(k) = C_factor * integral(func_exact, g_i, g_f);
% 
%     % --- B. 标准近似 (Standard) ---
%     % 假设 g / sqrt(g^2-1) ≈ 1
%     % 对应幂次 n = 0
%     val_std_f = get_R22_int(g_f, g_f, 0, Omega, K_sin);
%     val_std_i = get_R22_int(g_i, g_f, 0, Omega, K_sin);
%     I_standard(k) = C_factor * (val_std_f - val_std_i);
% 
%     % --- C. 修正近似 (Corrected) ---
%     % 假设 g / sqrt(g^2-1) ≈ 1 + 1/(2*g^2)
%     % 第一项是 Standard (n=0)
%     % 第二项对应幂次 n = -2
%     val_cor_f = get_R22_int(g_f, g_f, -2, Omega, K_sin);
%     val_cor_i = get_R22_int(g_i, g_f, -2, Omega, K_sin);
% 
%     term_correction = 0.5 * (val_cor_f - val_cor_i);
%     I_corrected(k) = I_standard(k) + C_factor * term_correction;
% 
%     % 计算误差
%     err_std(k) = abs(I_standard(k) - I_exact(k))/abs(I_exact(k)) * 100;
%     err_cor(k) = abs(I_corrected(k) - I_exact(k))/abs(I_exact(k)) * 100;
% end
% 
% % ==========================================
% % 4. 绘图结果
% % ==========================================
% figure('Color', 'w', 'Position', [100, 100, 1000, 400]);
% 
% subplot(1, 2, 1);
% plot(gamma_in_list, I_exact, 'k-', 'LineWidth', 2); hold on;
% plot(gamma_in_list, I_standard, 'b--', 'LineWidth', 1.5);
% plot(gamma_in_list, I_corrected, 'r-.', 'LineWidth', 1.5);
% title('\int R_{22}^{s->s_f} ds');
% xlabel('\gamma_{in}'); ylabel('Integral Value');
% legend('Exact', 'Standard (1)', 'Corrected (1 + 1/2\gamma^2)');
% grid on;
% 
% subplot(1, 2, 2);
% semilogy(gamma_in_list, err_std, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4); hold on;
% semilogy(gamma_in_list, err_cor, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
% yline(1, 'k:', '1% Threshold');
% title('relative err (%)');
% xlabel('\gamma_{in}');
% legend('Standard Error', 'Corrected Error');
% grid on;
% 
% fprintf('Gamma=2 时, 标准误差: %.2f%%, 修正误差: %.2f%%\n', err_std(1), err_cor(1));