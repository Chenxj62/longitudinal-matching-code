% 定义参数
%close all;
clear all;
k0=8.987551787e9;
e = 1.602176634e-19; % 电荷
R =10; % 弯曲半径
L_m = R*pi/9; % 最大长度
%% 
sigma = 3e-4; % 高斯分布的标准差
Ne=100e-12/e;

% 定义z/sigma的范围
z_sigma_range = linspace(-5, 5, 500); % 减少点数以加快计算速度
z_values = z_sigma_range * sigma;

% 初始化积分结果数组
integral_results = zeros(size(z_sigma_range));

% 对每个z值，计算dE/ds沿x从1到1.3的积分
for i = 1:length(z_sigma_range)
    z_fixed = z_values(i);
    
    % 定义x的范围，从1到1.3
    x_range = linspace(.1, 1.4, 50); % 积分用的点数
    dx=x_range(2)-x_range(1);
    
    % 初始化该z值下的积分结果
    integral_sum = 0;
    
    % 对每个x值计算dE/ds
    for j = 1:length(x_range)
        x = x_range(j);
        
        % 第一项: 4*lambda[s-Delta_s(L_m)]/(L_m+2x)
        Delta_s_L_m = (L_m^3 / (24 * R^2)) * ((L_m + 4*x) / (L_m + x));
        lambda_L_m = (1 / (sqrt(2*pi) * sigma)) * exp(-((z_fixed - Delta_s_L_m)^2) / (2*sigma^2));
        term1 = 4 * lambda_L_m / (L_m + 2*x);
        
        % 第二项: 积分项
        integration_points = 200; % 减少积分点数以加快计算
        l_values = linspace(L_m,0, integration_points);
        dl = L_m / (integration_points - 1);
        
        l_integral_sum = 0;
        for k = 1:integration_points
            l = l_values(k);
            
            % 避免l=0时可能的除零问题
            if l == 0
                integrand = 0;
            else
                Delta_s_l = (l^3 / (24 * R^2)) * ((l + 4*x) / (l + x));
                
                % 计算Delta_s(l)对l的偏导数
dDelta_s_dl = (l^2 / (8 * R^2)) * ((l + 2*x) / (l+x))^2;
                
                % 计算λ'[s-Δs(l)]
                z_shifted = z_fixed - Delta_s_l;
                lambda_derivative = -(z_shifted / (sqrt(2*pi) * sigma^3)) * exp(-(z_shifted^2) / (2*sigma^2));
                
                integrand = (4 / (l + 2*x)) * lambda_derivative * dDelta_s_dl;
            end
            
            % 梯形法则积分
            if k == 1 || k == integration_points
                weight = 0.5;
            else
                weight = 1.0;
            end
            
            l_integral_sum = l_integral_sum + 1 * integrand * dl;
        end
        
        % 计算该x值下的dE/ds
        dE_ds_value =  (term1 -l_integral_sum);
        
        % 使用矩形法则累加到x积分中
        if j == 1 || j == length(x_range)
            weight = 0.5;
        else
            weight = 1.0;
        end
        
        integral_sum = integral_sum +1 * dE_ds_value * dx;
    end
    
    % 存储该z值下的积分结果
    integral_results(i) = integral_sum;
end

% 绘制积分结果随z/sigma的变化
figure;
plot(z_sigma_range, Ne * k0 * e^2 *integral_results/1e9/e, 'LineWidth', 2);
xlabel('z/\sigma');
ylabel('∫(d\mathcal{E}/ds)dx [GeV·m/m]');
title('尾场沿x从1到1.3的积分随z/\sigma的变化');
grid on;