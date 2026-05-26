clear all;
close all;
x=[6.14809	-5.04067	1.30289	2.27912	0.2	4.44067	2.409	0.822748];
% ==== 参数定义 ====
sigz = 16e-6;      % 束长 σ_z
rho =4;           % 弯铁曲率半径 (m)
Lq = 0.2;          % 四极磁铁长度 (m)
allang = 15 * pi / 180; % 总偏转角 (弧度)
k = x(end);             % 优化变量 x(8)，用户可以调整
theta_final1 = allang * k / (2 + k); % 第一弯铁偏转角
theta_final2 = allang * 1 / (2 + k); % 第二弯铁偏转角
L_bend1 = rho * theta_final1;        % 第一弯铁长度
L_bend2 = rho * theta_final2;        % 第二弯铁长度
L_matching1 = 1;    % 第一匹配段长度
L_matching2 = 1;    % 第二匹配段长度
n_points = 2000;    % 计算点数
h = 0;              % 优化变量 h（可调整）



% ==== 匹配段矩阵 ====
% 第一匹配段矩阵
Q1 = Q(x(1), Lq);
Q2 = Q(x(2), Lq);
Q3 = Q(x(3), Lq);
D1 = D(x(4));
D2 = D(x(5));
D3 = D(x(6));
D4 = D(x(7));
M_matching1 = D4 * Q3 * D3 * Q2 * D2 * Q1 * D1;

% 第二匹配段矩阵
M_matching2 = D1 * Q1 * D2 * Q2 * D3 * Q3 * D4;

% ==== 创建完整束线的 s 数组 ====
L_total = 2 * L_bend2 + L_bend1 + L_matching1 + L_matching2;
s = linspace(0, L_total, n_points);

% 初始化积分和 R56 数组
R51_integral = zeros(size(s));
R52_integral = zeros(size(s));
R56_values = zeros(size(s));

% ==== 计算积分 ====
for i = 1:length(s)
    current_s = s(i);

    if current_s <= L_bend2 % 第一段弯铁
        theta_s = current_s / rho;
        M_total = bend_matrix(theta_s, rho);
        R56_values(i) = M_total(5, 6);

        if i > 1
            ds = s(i) - s(i - 1);
            C = 1 / (1 + h * R56_values(i));
            C_power = C^(4 / 3);
            R51_integral(i) = R51_integral(i - 1) + M_total(5, 1) * C_power * ds;
            R52_integral(i) = R52_integral(i - 1) + M_total(5, 2) * C_power * ds;
        end

    elseif current_s <= L_bend2 + L_matching1 % 第一匹配段
        idx_end_bend1 = find(s <= L_bend2, 1, 'last');
        R51_integral(i) = R51_integral(idx_end_bend1);
        R52_integral(i) = R52_integral(idx_end_bend1);
        M_total = M_matching1 * bend_matrix(theta_final2, rho);
        R56_values(i) = M_total(5, 6);

    elseif current_s <= L_bend1 + L_bend2 + L_matching1 % 第二段弯铁
        local_s = current_s - (L_bend2 + L_matching1);
        theta_s = local_s / rho;
        M_before = M_matching1 * bend_matrix(theta_final2, rho);
        M_current = bend_matrix(theta_s, rho);
        M_total = M_current * M_before;
        R56_values(i) = M_total(5, 6);

        idx_end_matching1 = find(s <= L_bend2 + L_matching1, 1, 'last');

        if i > idx_end_matching1
            ds = s(i) - s(i - 1);
            C = 1 / (1 + h * R56_values(i));
            C_power = C^(4 / 3);
            R51_integral(i) = R51_integral(i - 1) + M_total(5, 1) * C_power * ds;
            R52_integral(i) = R52_integral(i - 1) + M_total(5, 2) * C_power * ds;
        else
            R51_integral(i) = R51_integral(idx_end_matching1);
            R52_integral(i) = R52_integral(idx_end_matching1);
        end

    elseif current_s <= L_bend1 + L_bend2 + L_matching1 + L_matching2 % 第二匹配段
        idx_end_bend2 = find(s <= L_bend1 + L_bend2 + L_matching1, 1, 'last');
        R51_integral(i) = R51_integral(idx_end_bend2);
        R52_integral(i) = R52_integral(idx_end_bend2);
        M_total = M_matching2 * bend_matrix(theta_final1, rho) * M_matching1 * bend_matrix(theta_final2, rho);
        R56_values(i) = M_total(5, 6);

   else % 第三段弯铁
    local_s = current_s - (L_bend1 + L_bend2 + L_matching1 + L_matching2);
    theta_s = local_s / rho;
    M_before = M_matching2 * bend_matrix(theta_final1, rho) * M_matching1 * bend_matrix(theta_final2, rho);
    M_current = bend_matrix(theta_s, rho);
    M_total = M_current * M_before;
    R56_values(i) = M_total(5, 6);

    % 取前一个点的积分值
    idx_end_matching2 = find(s <= L_bend1 + L_bend2 + L_matching1 + L_matching2, 1, 'last');
    
    ds = s(i) - s(i - 1);
    C = 1 / (1 + h * R56_values(i));
    C_power = C^(4 / 3);
    R51_integral(i) = R51_integral(i - 1) + M_total(5, 1) * C_power * ds;
    R52_integral(i) = R52_integral(i - 1) + M_total(5, 2) * C_power * ds;

    end
end

% ==== 绘制积分图 ====
figure;
plot(s, R51_integral, 'r', 'LineWidth', 1.5, 'DisplayName', 'R51 Integral');
hold on;
plot(s, R52_integral, 'b', 'LineWidth', 1.5, 'DisplayName', 'R52 Integral');
xlabel('s (m)');
ylabel('Integral Value');
legend('show');
grid on;
title('R51 and R52 Integrals');

% 第一段弯铁传输矩阵
M_bend1 = bend_matrix(theta_final1, rho);

% 第二段弯铁传输矩阵
M_bend2 = bend_matrix(theta_final2, rho);

% 总传输矩阵
M_total = M_bend2 * M_matching2 * M_bend1 * M_matching1 * M_bend2;

% ==== 提取 tot56 和 tot16 ====
tot56 = M_total(5, 6); % 总传输矩阵中第 5 行第 6 列
tot16 = M_total(1, 6); % 总传输矩阵中第 1 行第 6 列
tot11 = M_total(1, 1) % 总传输矩阵中第 5 行第 6 列
tot12 = M_total(1, 2); % 总传输矩阵中第 1 行第 6 列

% ==== 输出结果 ====
fprintf('Total M56 (tot56): %.6e\n', tot56);
fprintf('Total M16 (tot16): %.6e\n', tot16);
fprintf('period beta: %.6e\n', sqrt(tot12^2/(1-tot11^2)));
function M = bend_matrix(theta, rho)
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    M = [cos_theta, rho * sin_theta, 0, 0, 0, rho * (1 - cos_theta);
         -sin_theta / rho, cos_theta, 0, 0, 0, sin_theta;
         0, 0, 1, rho * theta, 0, 0;
         0, 0, 0, 1, 0, 0;
         -sin_theta, -rho * (1 - cos_theta), 0, 0, 1, -rho * (theta - sin_theta);
         0, 0, 0, 0, 0, 1];
end