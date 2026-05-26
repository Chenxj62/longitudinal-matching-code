function [final_beta, evolution, periodic_twiss] = calcBetaEvolution(beamline, initial_energy, initial_beta, plot_flag, calc_periodic)
% calcBetaEvolution - 计算束线Beta函数演化和周期解（改进版，含R16跟踪）
%
% 输入参数:
%   beamline - 元件结构体数组
%   initial_energy - 初始能量(MeV)
%   initial_beta - 初始Beta参数 [betax, betay, alphax, alphay]
%   plot_flag - 是否绘图 (可选，默认false)
%   calc_periodic - 是否计算周期解 (可选，默认true)
%
% 输出参数:
%   final_beta - 最终Beta参数 [betax, betay, alphax, alphay]
%   evolution - 束线Beta演化记录结构体
%   periodic_twiss - 周期Twiss参数 [betax, betay, alphax, alphay]

%% 参数处理
if nargin < 4
    plot_flag = false;
end
if nargin < 5
    calc_periodic = false;
end

%% 先计算周期解（如果需要）
periodic_twiss = [];
use_periodic_initial = false;

if calc_periodic
    % 快速计算总传输矩阵以获得周期解
    [M_total_x, M_total_y] = calculateTotalMatrix(beamline);
    
    [beta_x_per, beta_y_per, alpha_x_per, alpha_y_per] = ...
        calculatePeriodicTwiss(M_total_x, M_total_y);
    
    has_x_solution = ~isempty(beta_x_per);
    has_y_solution = ~isempty(beta_y_per);
    
    if has_x_solution || has_y_solution
        if has_x_solution && has_y_solution
            periodic_twiss = [beta_x_per, beta_y_per, alpha_x_per, alpha_y_per];
            use_periodic_initial = true;
        elseif has_x_solution
            periodic_twiss = [beta_x_per, NaN, alpha_x_per, NaN];
            % 只有x方向有周期解，使用x方向周期解 + y方向初值
            initial_beta = [beta_x_per, initial_beta(2), alpha_x_per, initial_beta(4)];
            use_periodic_initial = true;
        else
            periodic_twiss = [NaN, beta_y_per, NaN, alpha_y_per];
            % 只有y方向有周期解，使用y方向周期解 + x方向初值
            initial_beta = [initial_beta(1), beta_y_per, initial_beta(3), alpha_y_per];
            use_periodic_initial = true;
        end
        
        % 如果两个方向都有周期解，使用周期解作为初值
        if has_x_solution && has_y_solution
            initial_beta = periodic_twiss;
        end
        
        % 输出信息
        trace_x = trace(M_total_x);
        trace_y = trace(M_total_y);
        
        fprintf('\n=== 周期/对称解计算结果 ===\n');
        if has_x_solution
            type_x = getStabilityType(trace_x);
            fprintf('X方向: β=%.4f m, α=%.4f (%s, trace=%.4f)\n', ...
                beta_x_per, alpha_x_per, type_x, trace_x);
        else
            fprintf('X方向: 无解 (trace=%.4f)\n', trace_x);
        end
        
        if has_y_solution
            type_y = getStabilityType(trace_y);
            fprintf('Y方向: β=%.4f m, α=%.4f (%s, trace=%.4f)\n', ...
                beta_y_per, alpha_y_per, type_y, trace_y);
        else
            fprintf('Y方向: 无解 (trace=%.4f)\n', trace_y);
        end
        
        if use_periodic_initial
            fprintf('\n>>> 使用周期解作为初始条件进行演化计算\n');
        end
        fprintf('===========================\n\n');
    else
        warning('无法计算有效的Twiss参数，使用输入初值');
    end
end

%% 初始化演化记录
num_elements = length(beamline);
num_slices = 10;  % 每个元件的切片数
max_points = num_elements * (num_slices + 2) + 100;

evolution = struct();
evolution.s = zeros(max_points, 1);
evolution.energy = zeros(max_points, 1);
evolution.beta_x = zeros(max_points, 1);
evolution.beta_y = zeros(max_points, 1);
evolution.alpha_x = zeros(max_points, 1);
evolution.alpha_y = zeros(max_points, 1);
evolution.gamma_x = zeros(max_points, 1);
evolution.gamma_y = zeros(max_points, 1);
evolution.element = cell(max_points, 1);

% 只在绘图时初始化R16
if plot_flag
    evolution.R16 = zeros(max_points, 1);
end

%% 初始值设置
evolution.s(1) = 0;
evolution.energy(1) = initial_energy;
evolution.beta_x(1) = initial_beta(1);
evolution.beta_y(1) = initial_beta(2);
evolution.alpha_x(1) = initial_beta(3);
evolution.alpha_y(1) = initial_beta(4);
evolution.gamma_x(1) = (1 + initial_beta(3)^2) / initial_beta(1);
evolution.gamma_y(1) = (1 + initial_beta(4)^2) / initial_beta(2);
evolution.element{1} = 'Start';

if plot_flag
    evolution.R16(1) = 0;
end

%% 当前状态变量
current_beta_x = initial_beta(1);
current_beta_y = initial_beta(2);
current_alpha_x = initial_beta(3);
current_alpha_y = initial_beta(4);
current_energy = initial_energy;
current_s = 0;
point_index = 2;

%% 总传输矩阵（用于验证）
M_total_x = eye(2);
M_total_y = eye(2);

% 只在绘图时累积完整6D矩阵
if plot_flag
    R_accumulated = eye(6);
end

%% 遍历束线元件
for i = 1:length(beamline)
    element = beamline(i);
    
    % 字段检查
    if ~isfield(element, 'type')
        warning('元素 #%d 没有type字段，跳过', i);
        continue;
    end
    
    if ~isfield(element, 'length')
        element.length = 0;
    end
    
    if ~isfield(element, 'name')
        element.name = sprintf('Unnamed_%s_%d', element.type, i);
    end
    
    %% 元件处理 - 统一逻辑
    if element.length <= 0
        % 零长度元件 - 直接使用传输矩阵
        R = getElementMatrix(element);
        M_x = R(1:2, 1:2);
        M_y = R(3:4, 3:4);
        
        % 只在绘图时累积完整矩阵
        if plot_flag
            R_accumulated = R * R_accumulated;
        end
        
        % 更新Beta参数
        [current_beta_x, current_alpha_x] = updateTwissParameters(...
            current_beta_x, current_alpha_x, M_x);
        [current_beta_y, current_alpha_y] = updateTwissParameters(...
            current_beta_y, current_alpha_y, M_y);
        
        % 更新能量
        if abs(R(6,6) - 1) > 1e-10
            current_energy = current_energy / R(6,6);
        end
        
        % 累积总传输矩阵
        M_total_x = M_x * M_total_x;
        M_total_y = M_y * M_total_y;
        
        % 记录终点（与起点相同位置）
        evolution.s(point_index) = current_s;
        evolution.energy(point_index) = current_energy;
        evolution.beta_x(point_index) = current_beta_x;
        evolution.beta_y(point_index) = current_beta_y;
        evolution.alpha_x(point_index) = current_alpha_x;
        evolution.alpha_y(point_index) = current_alpha_y;
        evolution.gamma_x(point_index) = (1 + current_alpha_x^2) / current_beta_x;
        evolution.gamma_y(point_index) = (1 + current_alpha_y^2) / current_beta_y;
        evolution.element{point_index} = element.name;
        
        if plot_flag
            evolution.R16(point_index) = R_accumulated(1,6);
        end
        
        point_index = point_index + 1;
        
    else
        % 有长度元件 - 使用切片法
        slice_length = element.length / num_slices;
        
        % 切片状态变量
        slice_beta_x = current_beta_x;
        slice_beta_y = current_beta_y;
        slice_alpha_x = current_alpha_x;
        slice_alpha_y = current_alpha_y;
        slice_energy = current_energy;
        
        % 元件累积传输矩阵
        M_element_x = eye(2);
        M_element_y = eye(2);
        
        % 处理每个切片
        for slice = 1:num_slices
            % 创建切片元件
            slice_element = element;
            slice_element.length = slice_length;
            
            % 获取切片传输矩阵
            R_slice = getElementMatrix(slice_element);
            M_x = R_slice(1:2, 1:2);
            M_y = R_slice(3:4, 3:4);
            
            % 只在绘图时累积完整矩阵
            if plot_flag
                R_accumulated = R_slice * R_accumulated;
            end
            
            % 更新Beta参数
            [slice_beta_x, slice_alpha_x] = updateTwissParameters(...
                slice_beta_x, slice_alpha_x, M_x);
            [slice_beta_y, slice_alpha_y] = updateTwissParameters(...
                slice_beta_y, slice_alpha_y, M_y);
            
            % 更新能量
            if abs(R_slice(6,6) - 1) > 1e-10
                slice_energy = slice_energy / R_slice(6,6);
            end
            
            % 累积元件传输矩阵
            M_element_x = M_x * M_element_x;
            M_element_y = M_y * M_element_y;
            
            % 记录所有切片点（包括最后一个）
            evolution.s(point_index) = current_s + slice * slice_length;
            evolution.energy(point_index) = slice_energy;
            evolution.beta_x(point_index) = slice_beta_x;
            evolution.beta_y(point_index) = slice_beta_y;
            evolution.alpha_x(point_index) = slice_alpha_x;
            evolution.alpha_y(point_index) = slice_alpha_y;
            evolution.gamma_x(point_index) = (1 + slice_alpha_x^2) / slice_beta_x;
            evolution.gamma_y(point_index) = (1 + slice_alpha_y^2) / slice_beta_y;
            evolution.element{point_index} = sprintf('%s_slice%d', element.name, slice);
            
            if plot_flag
                evolution.R16(point_index) = R_accumulated(1,6);
            end
            
            point_index = point_index + 1;
        end
        
        % 累积总传输矩阵
        M_total_x = M_element_x * M_total_x;
        M_total_y = M_element_y * M_total_y;
        
        % 更新当前状态
        current_beta_x = slice_beta_x;
        current_beta_y = slice_beta_y;
        current_alpha_x = slice_alpha_x;
        current_alpha_y = slice_alpha_y;
        current_energy = slice_energy;
        current_s = current_s + element.length;
    end
end

%% 裁剪数组
evolution.s = evolution.s(1:point_index-1);
evolution.energy = evolution.energy(1:point_index-1);
evolution.beta_x = evolution.beta_x(1:point_index-1);
evolution.beta_y = evolution.beta_y(1:point_index-1);
evolution.alpha_x = evolution.alpha_x(1:point_index-1);
evolution.alpha_y = evolution.alpha_y(1:point_index-1);
evolution.gamma_x = evolution.gamma_x(1:point_index-1);
evolution.gamma_y = evolution.gamma_y(1:point_index-1);
evolution.element = evolution.element(1:point_index-1);

if plot_flag
    evolution.R16 = evolution.R16(1:point_index-1);
end

%% 返回最终Beta参数
final_beta = [current_beta_x, current_beta_y, current_alpha_x, current_alpha_y];

%% 绘图
if plot_flag
    plotBetaEvolution(evolution, periodic_twiss, use_periodic_initial);
end

end

%% ==================== 辅助函数 ====================

function [M_total_x, M_total_y] = calculateTotalMatrix(beamline)
% 快速计算总传输矩阵（不记录演化）

M_total_x = eye(2);
M_total_y = eye(2);

for i = 1:length(beamline)
    element = beamline(i);
    
    if ~isfield(element, 'type')
        continue;
    end
    
    R = getElementMatrix(element);
    M_x = R(1:2, 1:2);
    M_y = R(3:4, 3:4);
    
    M_total_x = M_x * M_total_x;
    M_total_y = M_y * M_total_y;
end

end

function [beta_new, alpha_new] = updateTwissParameters(beta, alpha, M)
% 使用传输矩阵更新Twiss参数
% 输入:
%   beta, alpha - 当前Twiss参数
%   M - 2x2传输矩阵
% 输出:
%   beta_new, alpha_new - 更新后的Twiss参数

% 计算gamma
gamma = (1 + alpha^2) / beta;

% 构建Twiss矩阵
T = [beta, -alpha; -alpha, gamma];

% 应用传输: T_new = M * T * M'
T_new = M * T * M';

% 提取新参数
beta_new = T_new(1,1);
alpha_new = -T_new(1,2);

% 验证正定性
if beta_new <= 0
    warning('Beta函数变为非正值: %.6e，可能存在计算错误', beta_new);
    beta_new = abs(beta_new); % 取绝对值作为应急处理
end

end

function type_str = getStabilityType(trace_val)
% 根据trace值判断稳定性类型
if abs(trace_val) < 2
    type_str = '稳定-周期解';
elseif abs(trace_val) > 2
    type_str = '不稳定-对称解';
else
    type_str = '临界';
end
end

function plotBetaEvolution(evolution, periodic_twiss, is_periodic_evolution)
% 绘制Beta函数和R16演化图（双y轴）
% is_periodic_evolution: 标识当前绘制的是否为周期解演化

% 颜色定义
deep_blue = [0, 0, 0.7];
deep_red = [0.7, 0, 0];
light_blue = [0.3, 0.3, 1];
light_red = [1, 0.3, 0.3];
black = [0, 0, 0];

% 创建图形
figure('Position', [100, 100, 1200, 600], 'Color', 'w');
set(0, 'DefaultAxesFontName', 'Cambria');
set(0, 'DefaultTextFontName', 'Cambria');

%% 子图1: Beta函数和R16（双y轴）
subplot(2,1,1);

% 创建左轴（Beta）
yyaxis left
hold on;

% 绘制Beta演化
h1 = plot(evolution.s, evolution.beta_x, 'Color', deep_blue, 'LineWidth', 2);
h2 = plot(evolution.s, evolution.beta_y, 'Color', deep_red, 'LineWidth', 2,'LineStyle','-');

% 如果有周期解且当前不是周期解演化，绘制参考线
if ~isempty(periodic_twiss) && ~is_periodic_evolution
    if ~isnan(periodic_twiss(1))
        yline(periodic_twiss(1), '-', 'Color', light_blue, 'LineWidth', 1.5, ...
            'Label', sprintf('β_x^{per}=%.2f', periodic_twiss(1)));
    end
    if ~isnan(periodic_twiss(2))
        yline(periodic_twiss(2), '-', 'Color', light_red, 'LineWidth', 1.5, ...
            'Label', sprintf('β_y^{per}=%.2f', periodic_twiss(2)));
    end
end

ylabel('Beta (m)', 'FontSize', 12);
ylim([0, max([max(evolution.beta_x), max(evolution.beta_y)]) * 1.1]);
ax = gca;
ax.YColor = [0, 0, 0];

% 创建右轴（R16）
yyaxis right
hold on;

% 绘制R16演化
h3 = plot(evolution.s, evolution.R16, 'Color', black, 'LineWidth', 2, 'LineStyle', '-');

ylabel('R_{16} (m)', 'FontSize', 12);
ax = gca;
ax.YColor = black;

% 设置x轴
xlabel('s (m)', 'FontSize', 12);
xlim([0, max(evolution.s)]);

% 标题和图例
if is_periodic_evolution
    title('Beta Functions and Dispersion R_{16} (Periodic Solution Evolution)', ...
        'FontSize', 13, 'FontWeight', 'bold');
else
    title('Beta Functions and Dispersion R_{16}', 'FontSize', 13, 'FontWeight', 'bold');
end

legend([h1, h2, h3], {'\beta_x', '\beta_y', 'R_{16}'}, ...
    'Location', 'best', 'FontSize', 10);
grid on;
box on;

%% 子图2: Alpha函数
subplot(2,1,2);
hold on;

% 绘制Alpha演化
plot(evolution.s, evolution.alpha_x, 'Color', deep_blue, 'LineWidth', 2);
plot(evolution.s, evolution.alpha_y, 'Color', deep_red, 'LineWidth', 2,'LineStyle','-');

% 如果有周期解且当前不是周期解演化，绘制参考线
if ~isempty(periodic_twiss) && ~is_periodic_evolution
    if ~isnan(periodic_twiss(3))
        yline(periodic_twiss(3), '--', 'Color', light_blue, 'LineWidth', 1.5, ...
            'Label', sprintf('α_x^{per}=%.2f', periodic_twiss(3)));
    end
    if ~isnan(periodic_twiss(4))
        yline(periodic_twiss(4), '--', 'Color', light_red, 'LineWidth', 1.5, ...
            'Label', sprintf('α_y^{per}=%.2f', periodic_twiss(4)));
    end
end

% 添加零线参考
yline(0, 'k:', 'LineWidth', 0.5);

xlabel('s (m)', 'FontSize', 12);
ylabel('Alpha', 'FontSize', 12);

if is_periodic_evolution
    title('Alpha Functions (Periodic Solution Evolution)', 'FontSize', 13, 'FontWeight', 'bold');
else
    title('Alpha Functions', 'FontSize', 13, 'FontWeight', 'bold');
end

legend('\alpha_x', '\alpha_y', 'Location', 'best', 'FontSize', 10);
grid on;
xlim([0, max(evolution.s)]);
box on;

end

%% ==================== 周期解计算函数 ====================

function [beta_x, beta_y, alpha_x, alpha_y] = calculatePeriodicTwiss(M_x, M_y)
% 计算周期/对称Twiss参数
[beta_x, alpha_x] = calculateSinglePlaneTwiss(M_x);
[beta_y, alpha_y] = calculateSinglePlaneTwiss(M_y);
end

function [beta, alpha] = calculateSinglePlaneTwiss(M)
% 单平面Twiss参数计算

beta = [];
alpha = [];

trace_val = trace(M);

%% 稳定系统 (|trace| < 2) - 解析解
if abs(trace_val) < 2
    try
        mu = acos(0.5 * trace_val);
        sin_mu = sin(mu);
        
        if abs(sin_mu) > 1e-10
            % 提取矩阵元素
            s = M(1,2);
            c = M(1,1);
            sp = M(2,2);
            
            % 计算beta和alpha
            beta = abs(s / sin_mu);
            alpha = (c - sp) / (2 * sin_mu);
            
            % 验证
            if ~verifyTwissSolution(beta, alpha, M, 'periodic')
                beta = [];
                alpha = [];
            end
        end
    catch
        % 计算失败
    end
    
%% 不稳定系统 (|trace| > 2) - 数值解
elseif abs(trace_val) > 2
    % 方法1: 解析方法（基于特征值）
    [beta, alpha] = solveSymmetricAnalytic(M);
    
    % 如果解析方法失败，尝试数值方法
    if isempty(beta) && exist('fsolve', 'file')
        [beta, alpha] = solveSymmetricNumerical(M);
    end
end

end

function [beta, alpha] = solveSymmetricAnalytic(M)
% 解析方法求解对称解

beta = [];
alpha = [];

try
    % 对于对称解: M*T*M' = T_flip，其中T_flip = [β, α; α, γ]
    % 这导致方程: M12*α^2 + (M11-M22)*α - M21 = 0
    
    A = M(1,2);
    B = M(1,1) - M(2,2);
    C = -M(2,1);
    
    if abs(A) < 1e-10
        return;
    end
    
    discriminant = B^2 - 4*A*C;
    
    if discriminant >= 0
        alpha1 = (-B + sqrt(discriminant))/(2*A);
        alpha2 = (-B - sqrt(discriminant))/(2*A);
        
        % 计算beta
        if abs(alpha2 - alpha1) > 1e-10
            beta_val = abs(M(1,2) / (alpha2 - alpha1));
            
            % 验证两个alpha值，选择有效的
            if verifyTwissSolution(beta_val, alpha1, M, 'symmetric')
                beta = beta_val;
                alpha = alpha1;
            elseif verifyTwissSolution(beta_val, alpha2, M, 'symmetric')
                beta = beta_val;
                alpha = alpha2;
            end
        end
    end
catch
    % 解析方法失败
end

end

function [beta, alpha] = solveSymmetricNumerical(M)
% 数值方法求解对称解

beta = [];
alpha = [];

try
    % 初始猜测
    s = M(1,2);
    if abs(s) > 1e-10
        beta0 = abs(s);
        alpha0 = 0;
    else
        beta0 = 1;
        alpha0 = 0;
    end
    
    % 定义方程
    eqn = @(x) symmetricTwissResidual(x, M);
    
    % 设置求解选项
    options = optimoptions('fsolve', ...
        'Display', 'off', ...
        'Algorithm', 'levenberg-marquardt', ...
        'MaxIterations', 10000, ...
        'MaxFunctionEvaluations', 10000, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-12);
    
    % 求解
    [result, ~, exitflag] = fsolve(eqn, [beta0, alpha0], options);
    
    if exitflag > 0 && result(1) > 0
        if verifyTwissSolution(result(1), result(2), M, 'symmetric')
            beta = result(1);
            alpha = result(2);
        end
    end
catch
    % 数值方法失败
end

end

function residual = symmetricTwissResidual(params, M)
% 对称解残差函数

beta = params(1);
alpha = params(2);

if beta <= 0
    residual = [1e10; 1e10];
    return;
end

gamma = (1 + alpha^2) / beta;
T = [beta, -alpha; -alpha, gamma];

T_transformed = M * T * M';

% 对称条件: T_transformed = [β, α; α, γ]
residual = [
    T_transformed(1,1) - beta;
    T_transformed(1,2) + alpha  % 注意符号
];

end

function is_valid = verifyTwissSolution(beta, alpha, M, solution_type)
% 验证Twiss参数解的有效性

is_valid = false;

if beta <= 0 || isnan(beta) || isnan(alpha)
    return;
end

try
    gamma = (1 + alpha^2) / beta;
    T = [beta, -alpha; -alpha, gamma];
    T_transformed = M * T * M';
    
    beta_t = T_transformed(1,1);
    alpha_t = -T_transformed(1,2);
    
    tolerance = 1e-6;
    
    switch solution_type
        case 'periodic'
            % 周期解: T_transformed = T
            is_valid = (abs(beta_t - beta) < tolerance) && ...
                       (abs(alpha_t - alpha) < tolerance);
        case 'symmetric'
            % 对称解: T_transformed = T_flip
            is_valid = (abs(beta_t - beta) < tolerance) && ...
                       (abs(alpha_t + alpha) < tolerance);
    end
catch
    is_valid = false;
end

end