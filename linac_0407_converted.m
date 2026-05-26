% Converted from linac_0407.lte (elegant lattice)
% Kickers, monitors, TDS, bend (no angle) merged into adjacent drifts

clear all;

close all;

set(0, 'DefaultAxesFontName', 'Cambria');
set(0, 'DefaultTextFontName', 'Cambria');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 10);
cons=[];

E_target =10.9e6;   % 目标能量 10 MeV (eV)



% ========== Match Section ==========
% Drifts
ELE.LM_D0_1   = createElement('drift', 0.005, 0, 'LM_D0_1');
ELE.LM_K01_1  = createElement('drift', 0.1, 0, 'LM_K01_1');     % ekicker
ELE.LM_K01_2  = createElement('drift', 0.1, 0, 'LM_K01_2');     % ekicker
ELE.LM_D0_2   = createElement('drift', 0.22, 0, 'LM_D0_2');
ELE.LM_YAG01  = createElement('drift', 0.05, 0, 'LM_YAG01');    % monitor
ELE.LM_D0_3   = createElement('drift', 0.05, 0, 'LM_D0_3');
ELE.LM_D1     = createElement('drift', 0.3, 0, 'LM_D1');
ELE.LM_D2     = createElement('drift', 0.3, 0, 'LM_D2');
ELE.LM_D3     = createElement('drift', 0.3, 0, 'LM_D3');
ELE.LM_D4     = createElement('drift', 0.1, 0, 'LM_D4');
ELE.LM_TDS    = createElement('drift', 0.5, 0, 'LM_TDS');       % TDS rfca, freq=2.856GHz
ELE.LM_D5_1   = createElement('drift', 0.1, 0, 'LM_D5_1');
ELE.LM_YAG02  = createElement('drift', 0.15, 0, 'LM_YAG02');    % monitor
ELE.LM_D5_2   = createElement('drift', 0.2, 0, 'LM_D5_2');
ELE.LM_B1     = createElement('drift', 0.22, 0, 'LM_B1');       % csbend, no angle
ELE.LM_D6_1   = createElement('drift', 1.07, 0, 'LM_D6_1');
ELE.LM_K03_1  = createElement('drift', 0.1, 0, 'LM_K03_1');     % ekicker
ELE.LM_K03_2  = createElement('drift', 0.1, 0, 'LM_K03_2');     % ekicker
ELE.LM_D6_2   = createElement('drift', 0.04, 0, 'LM_D6_2');
ELE.LM_D7     = createElement('drift', 0.3, 0, 'LM_D7');
ELE.LM_D8     = createElement('drift', 0.3, 0, 'LM_D8');
ELE.LM_D9_1   = createElement('drift', 0.1, 0, 'LM_D9_1');
ELE.LM_K04_1  = createElement('drift', 0.1, 0, 'LM_K04_1');     % ekicker
ELE.LM_K04_2  = createElement('drift', 0.1, 0, 'LM_K04_2');     % ekicker
ELE.LM_D9_2   = createElement('drift', 2.40316, 0, 'LM_D9_2');

% Quadrupoles (match section)
ELE.LM_Q1 = createElement('quadrupole', 0.12, -10.5119393, 'LM_Q1');
ELE.LM_Q2 = createElement('quadrupole', 0.12,   4.4322476, 'LM_Q2');
ELE.LM_Q3 = createElement('quadrupole', 0.12,   5.1127467, 'LM_Q3');
ELE.LM_Q4 = createElement('quadrupole', 0.12,   4.5482837, 'LM_Q4');
ELE.LM_Q5 = createElement('quadrupole', 0.12, -11.9836732, 'LM_Q5');
ELE.LM_Q6 = createElement('quadrupole', 0.12,   9.3179641, 'LM_Q6');
ELE.LM_Q7 = createElement('quadrupole', 0.12,   2.9055796, 'LM_Q7');

% ========== Linac Section ==========
% Drifts
ELE.DCAVMAP_A   = createElement('drift', 0.173128, 0, 'DCAVMAP_A');
ELE.DCAVMAP_B   = createElement('drift', 0.173128, 0, 'DCAVMAP_B');
ELE.DQA1        = createElement('drift', 0.115, 0, 'DQA1');
ELE.DQA2        = createElement('drift', 0.17, 0, 'DQA2');
ELE.DQB         = createElement('drift', 0.663, 0, 'DQB');

% Cavity (all 32 cavities identical)
% gradient = 15.625 MV / 1.03774 m = 15.058 MV/m, 9 cells, 1300 MHz
% phase = 97 deg (elegant: V*sin(phase), on-crest=90, i.e. 7 deg off-crest)
ELE.CAV = createElement('cavity', 0, [15.058, 9, 1300, 0], 'CAV');

% Quadrupoles (linac inter-cryomodule)
ELE.QCM01 = createElement('quadrupole', 0.2,  0.63672, 'QCM01');
ELE.QCM02 = createElement('quadrupole', 0.2, -0.39600, 'QCM02');
ELE.QCM03 = createElement('quadrupole', 0.2,  0.22418, 'QCM03');

% ========== Build beamline ======== Match section
MATCH_SEQ = [ELE.LM_D0_1, ELE.LM_K01_1, ELE.LM_K01_2, ...
    ELE.LM_D0_2, ELE.LM_YAG01, ELE.LM_D0_3, ...
    ELE.LM_Q1, ELE.LM_D1, ...
    ELE.LM_Q2, ELE.LM_D2, ...
    ELE.LM_Q3, ELE.LM_D3, ...
    ELE.LM_Q4, ELE.LM_D4, ...
    ELE.LM_TDS, ELE.LM_D5_1, ELE.LM_YAG02, ELE.LM_D5_2, ...
    ELE.LM_B1, ELE.LM_D6_1, ELE.LM_K03_1, ELE.LM_K03_2, ELE.LM_D6_2, ...
    ELE.LM_Q5, ELE.LM_D7, ...
    ELE.LM_Q6, ELE.LM_D8, ...
    ELE.LM_Q7, ELE.LM_D9_1, ELE.LM_K04_1, ELE.LM_K04_2, ELE.LM_D9_2];

% Cavity cluster: 8 cavities
CAV_CLUSTER = [ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B];

% Inter-cryomodule quad sections
QUAD_SEC_1 = [ELE.DQA1, ELE.DQA2, ELE.QCM01, ELE.DQB];
QUAD_SEC_2 = [ELE.DQA1, ELE.DQA2, ELE.QCM02, ELE.DQB];
QUAD_SEC_3 = [ELE.DQA1, ELE.DQA2, ELE.QCM03, ELE.DQB];
END_SEC    = [ELE.DQA1, ELE.DQA2];

% Linac section: 4 cryomodules + 3 quad sections
LINAC_SEQ = [CAV_CLUSTER, QUAD_SEC_1, ...
    CAV_CLUSTER, QUAD_SEC_2, ...
    CAV_CLUSTER, QUAD_SEC_3, ...
    CAV_CLUSTER, END_SEC];

% Full beamline
beamline = [MATCH_SEQ , ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B, ...
    ELE.DCAVMAP_A, ELE.CAV, ELE.DCAVMAP_B];
r=calculateTotalMatrix(beamline,10.72/0.511);
[results] = TSCenv(beamline, E_target, 'fvv12_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values, emitx_proj_array, emity_proj_array, x_centroid, y_centroid]=SCEnv(ELE.MATCH_D0, init_params, 'fvv9mer_slice_params.txt');
%[s_array, sigx_array, sigy_array, k1_values] = SCEnv(beamline, init_params);%, 'beam_11_slice_params.txt');

% 第一张图：x方向束流包络(sigma_x)演化
% 创建图形窗口
figure('Position', [100, 100,600, 500]);



% ========== 第一张图：束流尺寸演化 ==========
subplot(3, 1, 1);
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
legend([h1, h2], {'analysis x', 'analysis y', 'RF Track','Ocelot', 'Astra'}, ...
       'Location', 'best', 'FontSize', 11);

% 设置刻度字体
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% ========== 第二张图：发射度演化 ==========
subplot(3, 1, 2);
% 绘制归一化发射度


% 绘制数据点 (data1 - 圆圈)

h3 = plot(results.s_array, results.gamma_array.* results.peak_slice_emitx , ...
          'r-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,x} (Calculation)');
hold on;
h4 = plot(results.s_array, results.gamma_array.* results.peak_slice_emity , ...
          'b-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,y} (Calculation)');

xlabel('s (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\epsilon_n (mm·mrad)', 'FontSize', 14, 'FontWeight', 'bold');
title('Normalized Emittance Evolution', 'FontSize', 15, 'FontWeight', 'bold');

% 添加网格和图例
grid off;
box on;

legend([h3, h4], {'analysis x', 'analysis y'}, ...
       'Location', 'best', 'FontSize', 11);

% 设置刻度字体
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
subplot(3, 1, 3);
h3 = plot(results.s_array, results.sigx_array.^2 ./ results.emitx_proj_array , ...
          'r-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,x} (Calculation)');
hold on;
h4 = plot(results.s_array, results.sigy_array.^2 ./ results.emity_proj_array , ...
          'b-', 'LineWidth', 2, 'DisplayName', '\epsilon_{n,y} (Calculation)');

xlabel('s (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('\beta (m)', 'FontSize', 14, 'FontWeight', 'bold');
%title('Normalized Emittance Evolution', 'FontSize', 15, 'FontWeight', 'bold');

% 添加网格和图例
grid off;
box on;

legend([h3, h4], {'analysis x', 'analysis y'});

% 设置刻度字体
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
 figure;
% 
% % 第一张图：x方向束流包络(sigma_x)演化
% subplot(4, 1, 1);
% plot(results.s_array, results.sigy_array * 1e3, 'b-', 'LineWidth', 1.5);
% xlabel('s (m)', 'FontSize', 12);
% ylabel('\sigma_x (mm)', 'FontSize', 12);
% title('Y sig', 'FontSize', 14);
% grid on;
% 
% % 第二张图：x方向质心位置(centx)演化
subplot(4, 1, 2);
plot(results.s_array, results.x_centroid * 1e3, 'r-', 'LineWidth', 1.5);
xlabel('s (m)', 'FontSize', 12);
ylabel('x_{cent} (mm)', 'FontSize', 12);
title('Y cent', 'FontSize', 14);
