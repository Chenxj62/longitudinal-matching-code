function [y,cons]=linacopt(x)


cons=[];

emitx = 5e-7;      % 归一化发射度 x (m·rad)
emity = 5e-7;      % 归一化发射度 y (m·rad)
betax = 23;        % beta函数 x (m)
betay = 143;        % beta函数 y (m)
alphax = 72;        % alpha函数 x
alphay = -4.3;        % alpha函数 y
E_target = 10.8e6;   % 目标能量 10 MeV (eV)
Q=5e-11;

% 物理常数

% 计算gamma因子






% Drift 元件
ELE.MATCH_D0 = createElement('drift', 0.5, 0, 'MATCH_D0');

ELE.D= createElement('drift', 0.1, 0, 'MATCH_D7');

% Quadrupole 元件 (长度, k1, 名称)
% k1 值放大10倍
ELE.MATCH_Q1 = createElement('quadrupole', 0.12, x(1), 'MATCH_Q1');
ELE.MATCH_Q2 = createElement('quadrupole', 0.12, x(2), 'MATCH_Q2');
ELE.MATCH_Q3 = createElement('quadrupole', 0.12, x(3), 'MATCH_Q3');
ELE.MATCH_Q4 = createElement('quadrupole', 0.12, x(4), 'MATCH_Q4');


% ========== Linac Section (加速段) ==========
ELE.LINAC_D0 = createElement('drift', 0.4, 0, 'LINAC_D1');
ELE.LINAC_D1 = createElement('drift', 0.4, 0, 'LINAC_D1');


% Quadrupole 元件 (长度, k1, 名称)
% k1 值放大10倍，第三个重复第一个
ELE.LINAC_Q1 = createElement('quadrupole', 0.2, x(5), 'LINAC_Q1');
ELE.LINAC_Q2 = createElement('quadrupole', 0.2, x(6), 'LINAC_Q2');
ELE.LINAC_Q5 = createElement('quadrupole', 0.2, x(8), 'LINAC_Q1');
ELE.LINAC_Q4 = createElement('quadrupole', 0.2, x(9), 'LINAC_Q2');
ELE.LINAC_Q3 = createElement('quadrupole', 0.2,x(7), 'LINAC_Q2');
% Cavity 元件 (长度, [梯度MV/m, CELL数, 频率MHz, 相位deg], 名称)
% 梯度: 15.15 MV/m, 1.038 m 长度, 1300 MHz, -10 deg
ELE.CAV = createElement('cavity', 0.2, [12.5, 9, 1300, 0], 'CAV');

% ========== 构建完整 beamline ==========
% Match Section (7个四极铁 + drifts)
MATCH_SEQ = [ELE.MATCH_D0, ELE.MATCH_Q1,ELE.MATCH_D0, ELE.MATCH_Q2,ELE.MATCH_D0, ELE.MATCH_Q3,ELE.MATCH_D0, ELE.MATCH_Q4];

% Cavity Cluster 1 (8个 cavity)
CAV_CLUSTER_1 = [
    ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0, ...
     ELE.LINAC_D0, ELE.CAV, ELE.LINAC_D0];


% 完整 Linac Section
LINAC_SEQ = [
    CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q1, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q2, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q3, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q4, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q5, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q4, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q3, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q2, ELE.LINAC_D1, ...
     CAV_CLUSTER_1, ELE.LINAC_D1, ELE.LINAC_Q1, ELE.LINAC_D1, ...
     CAV_CLUSTER_1];

 %data4 = load('nory.line_dat');
% 完整 beamline
beamline = [MATCH_SEQ,LINAC_SEQ];
% 初始参数
init_params = [betax, betay, alphax, alphay, emitx, emity, 1e-3,Q, E_target, 0.001];


[results] = TSCenv( beamline, init_params, 'fvv11_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values, emitx_proj_array, emity_proj_array, x_centroid, y_centroid]=SCEnv(ELE.MATCH_D0, init_params, 'fvv9mer_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values] = SCEnv(beamline, init_params);%, 'beam_11_slice_params.txt');
betax=max(results.sigx_array.^2./results.emitx_proj_array);
betay=max(results.sigy_array.^2./results.emity_proj_array);

enx=2000*results.emitx_proj_array(end);
eny=2000*results.emity_proj_array(end);

enxy=max(enx,eny);

y=[(betax+betay)/2,enxy];


subplot(2, 1, 1);
% 绘制计算结果

h1 = plot(results.s_array, results.sigx_array * 1e3, 'r-', 'LineWidth', 2);
hold on;
h2 = plot(results.s_array, results.sigy_array * 1e3, 'b-', 'LineWidth', 2);


% 设置坐标轴和标签
%xlim([0, 20]);
xlabel('s (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\sigma (mm)', 'FontSize', 14, 'FontWeight', 'bold');
%title('Beam Size Evolution: \sigma_{x} and \sigma_{y}', 'FontSize', 15, 'FontWeight', 'bold');

% 添加网格和图例
grid off;
box on;

% 设置刻度字体
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% ========== 第二张图：发射度演化 ==========
subplot(2, 1, 2);
% 绘制归一化发射度
 plot(results.s_array, results.gamma_array .* results.emitx_proj_array * 1e6, ...
          'r-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,x} (Calculation)');
hold on;
 plot(results.s_array, results.gamma_array .* results.emity_proj_array * 1e6, ...
          'b-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,y} (Calculation)');

% 绘制数据点 (data1 - 圆圈)



end