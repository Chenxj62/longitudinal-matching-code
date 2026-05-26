%%%%%%%%%%%%%%%BMAD plot


clear all;
close all;

figure;

pnx=400;
pny=400;%%%%%%%%%%%%%%%%%%%%%图像分辨率
beam_source='bo'; %%%%%%rft (rftrack) ast(astra) bo=bmad old 13列 bn bmad 新格式 sdds=elegant, 如果选2, 记得charge_ele
energy=1000.24;%%%%%%%%能量 MeV
nbin=500;   %%%%%%%纵向切片
c_sig=0;  %%%%%%%截断束团边缘, [%], 不截断选择0
cc=0;           %%%%%%%%%矫正中心cc=1, cc不为1则不矫正
charge=5e-11;%%%C  BMAD file 不用管
beamratio=4;
plotimag='longmat';  %x, y, xy, xz.(相空间自定义色温图系数tempnx(y)), 
               %z, 1x/1y=slice_emit, 2x/2y=center in z, 2xc/2yc=center in xxp/yyp, 3x/3y=slicetwi,

               % 4x/4y=phase elip n_sli=min[round(np/1000),200],(纵向根据粒子数改变切片数nbin),
               % se=slice energy spread  mex=mean emit x
               % we=write sdds, 
               % wgd=write genesis dist file, 
               % wgb=write genesis beamfile(nslice定义切片数)
               % twm=twiss match(自定义目标twiss--target_*), 

                %exb3/exb=expand beam (切片+采样)
                %ewb=elegant write bmad .charge_ele, 抽样比例, beam_source!!!!!
                %skew=计算偏度. mma=mismatching in x 记得maxk, 截断系数(0~1,不截断为1).
                %mmz= mismatch in z
                %inform: 基本信息
                %longmat=longitudinal match：
                %phi_0:定义中心加速相位,0为峰值加速，EA：加速能量，R56，T566，U5666
                %fit：是否拟合相空间，1：拟合 else否
    
infile='fvv11.beam';
outputfile='C:\Users\chenxj\Desktop\erlarc\9.20\fvv11.beam';


%%%%% twiss match, outfiel=outputfile  使用fsolve, 最大次数10000
target_betx=41.9;
target_bety=75;%17.04;
target_alpx=18.3;
target_alpy =-6.6;

% target_betx = 86.95;
% target_bety = 34;
% target_alpx = 6.498;
% target_alpy = 0.35;

target_betx=86.95;
target_bety=104;%17.04;
target_alpx=2;%6.498;
target_alpy =0.35;

target_betx = 4.3;  % 【请指定】X方向目标 Beta (m)
target_alpx = 10;   % 【请指定】X方向目标 Alpha
target_bety = -3.2;  % 【请指定】Y方向目标 Beta (m)
target_alpy = .5;   % 【请指定】Y方向目标 Alpha




%%%%%%%%%expand beam; 
n_beam=10;
perturbationScale=.01;
expand_slice=200; %only exb3, 扩束切片数
expand_edge=0 ; %only exb3, 边缘矫正比例 单位：%
%%%%%%%%%%%%%%%%x, y, z相空间色温图
tempnx=2000;
tempny=2000;
%%%%%%%%%%%%%%%%%genesis beamfile






%%%%%%%%%%%longmat
phi_0=0;   %%%%%%%%%单位：°,0：峰值加速
EA=0;  %%%MeV
fit=11;
zp0=1; %%% 如果为0，令z=0处pz=0，否则不矫正
R56=0.000;
T566=0.;%70366;
U5666=0;
freq=1300; %%%%MHz


 R560=0.00;%.149;
 T5660=-0.;%11;%70366;
 U56660=00;
maxr=.2;

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch beam_source

    case 'ast'

        nc=charge; %%%%%bunch charge [C]

% 添加比例参数
sample_ratio =beamratio;  % 设置抽样比例，如0.1表示抽取10%的粒子

% 读取文本文件数据
data = dlmread(infile);

% 根据条件筛选行
filteredData = data(data(:, 10) > 0, :);

% 提取前六列数据
extractedData = filteredData(:, [1, 4, 2, 5, 3, 6]);

np=size(extractedData,1)-1;

% 第一列与第一行的差值
extractedData(:, 1) = extractedData(:, 1) + extractedData(1, 1);

% 第二列与第一行的差值
extractedData(:, 2) = (extractedData(:, 2) + extractedData(1, 2))/ extractedData(1, 6);

% 第四列与第一行的差值并除以第一行的第六列
extractedData(:, 3) = extractedData(:, 3) + extractedData(1, 3);

% 第五列与第一行的差值并除以第一行的第六列
extractedData(:, 4) = (extractedData(:, 4) + extractedData(1, 4)) / extractedData(1, 6);

% 计算第六列的平均值
extractedData(:, 6) = extractedData(:, 6) + extractedData(1, 6);
energy = mean(extractedData(:, 6));
fprintf('束团能量为 %f MeV\n', energy/1e6);



% 计算第六列与平均值的相对误差
relativeError = (extractedData(:, 6) - energy) ./ energy;

energy=energy/1e6;

% 更新第六列元素
extractedData(:, 6) = relativeError;

extractedData(1,:)=[];
extractedData(:,1)=extractedData(:,1)-mean(extractedData(:,1));
extractedData(:,3)=extractedData(:,3)-mean(extractedData(:,3));
extractedData(:,2)=extractedData(:,2)-mean(extractedData(:,2));
extractedData(:,4)=extractedData(:,4)-mean(extractedData(:,4));

maxsigz=max(extractedData(:,5));
extractedData(extractedData(:,5)<-3*maxsigz,:)=[];
npastra=size(extractedData,1);

% 按比例随机抽取粒子

    sampledData = extractedData;
    np_sampled = npastra;



% 或者分开写：
x  = sampledData(:, 1);
xp = sampledData(:, 2);
y  = sampledData(:, 3);
yp = sampledData(:, 4);
z  = sampledData(:, 5);
p  = sampledData(:, 6);


case 'sdds'
    % ===== 解析SDDS头部 =====
    fid = fopen(infile, 'r');
    if fid == -1
        error('无法打开文件: %s', infile);
    end
    
    column_names = {};
    column_count = 0;
    has_charge_param = false;
    
    % 读取头部,提取列名和参数
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        
        % 跳过空行和注释
        if isempty(line) || startsWith(line, '!'), continue; end
        
        % 检测Charge参数
        if contains(line, '&parameter') && contains(lower(line), 'charge')
            has_charge_param = true;
        end
        
        % 提取列名
        if contains(line, '&column')
            % 提取 name=xxx (处理可能的空格)
            tokens = regexp(line, 'name\s*=\s*(\w+)', 'tokens');
            if ~isempty(tokens)
                column_count = column_count + 1;
                column_names{column_count} = tokens{1}{1};
            end
        end
        
        % 检测数据开始
        if contains(line, '&data')
            break;
        end
    end
    
    % ===== 读取数据段 =====
    % 第一行: Charge值 (如果有参数定义)
    if has_charge_param
        charge_line = fgetl(fid);
        nc = str2double(strtrim(charge_line));
        if isnan(nc)
            warning('无法读取Charge参数,使用默认值');
            nc = charge;
        end
    else
        nc = charge;
    end
    
    % 第二行: 粒子数
    np_line = fgetl(fid);
    np = str2double(strtrim(np_line));
    
    if isnan(np) || np <= 0
        error('无法读取粒子数,请检查文件格式');
    end
    
    % 构建格式字符串 (所有列都是浮点数)
    format_str = repmat('%f', 1, column_count);
    
    % 读取数据矩阵
    data_matrix = textscan(fid, format_str, np, 'CollectOutput', false);
    fclose(fid);
    
    % ===== 映射到变量 =====
    for i = 1:column_count
        col_name = lower(column_names{i});
        fprintf('i=%d\n',i)
        switch col_name
            case 't'
                z = data_matrix{i};
            case 'x'
                x = data_matrix{i};
            case 'xp'
                xp = data_matrix{i};
            case 'y'
                y = data_matrix{i};
            case 'yp'
                yp = data_matrix{i};
            case 'p'
                p = data_matrix{i};
    
            otherwise
                % 其他列也可以存储
                warning('未知列: %s', column_names{i});
        end
    end
    
    





    z=-(z-mean(z))*299792458;
    energy=mean(p)*0.511;
    fprintf('energy=%.2f\n',energy);
    p=1*(p-mean(p))/mean(p);

   


    case 'rft'
     [x, xp, y, yp, z, p]=textread(infile,'%n%n%n%n%n%n','headerlines',0);


np=size(x,1);

Ene=mean(p);

energy=Ene;

fprintf('束团能量为:%f MeV\n',Ene);

%% 创建纵向切片
x=x/1000;
xp=xp/1000;
y=y/1000;
yp=yp/1000;
z=z-mean(z);
z=-1*z/1000;

p=1*(p-mean(p))/mean(p);
nc=charge;
    case 'bo'
%arc1_20w.dat';
%infile='beam_arc1out.dat'

[np]=textread(infile,'%f%*s%*s%*s',1,'headerlines',3);        
[nc]=textread(infile,'%f%*s%*s%*s',1,'headerlines',6);
[x xp y yp z p]=textread(infile,'%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f',np,'headerlines',9);
%[x xp y yp z p]=textread(infile,'%n%n%n%n%n%n',np,'headerlines',9);


    case 'bn'
fid = fopen(infile,'r');

q_total = NaN;
n_particle = NaN;

% 逐行读 header，直到遇到列名行 "#!"
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('文件未找到数据段（未发现 #! 列名行）');
    end

    if contains(line, '# charge')
        q_total = sscanf(line, '# charge%*[^=]=%f');
    elseif contains(line, '# n_particle')
        n_particle = sscanf(line, '# n_particle%*[^=]=%f');
    end

    if startsWith(strtrim(line), '#!')
        break; % 下一行开始就是数据
    end
end

% 读取数据行（18列，按需跳过）
fmt = ['%f %f %f %f %f %f %f ' ...                 % index x px y py z pz
       '%*f %*f %*f %*f %*f %*f %*f %*f ' ...      % p0c s time spinx spiny spinz phasex phasey（跳过）
       '%s %f %s'];                                 % state ix_ele location

data = textscan(fid, fmt, 'Delimiter',' ', 'MultipleDelimsAsOne',true);
fclose(fid);

index = data{1};
x     = data{2};
xp    = data{3};
y     = data{4};
yp    = data{5};
z     = data{6};
p     = data{7};     % pz
state = data{8};     % 'Alive'
ix_ele = data{9};
location = data{10}; % 'Downstream_End'

np = numel(x);

% q_total 是总电荷，不要再乘 np
% 如果你想要“每个粒子电荷”，用：
nc=q_total* np;

end



if cc==1
x=x-mean(x);
xp=xp-mean(xp);
y=y-mean(y);
yp=yp-mean(yp);
z=z-mean(z);
else
    x=x;
    y=y;
end

phase=[x,xp,y,yp,z,p];


if abs(max(phase(:,5)))-abs(min(phase(:,5)))>2*std(z)

rowsToDelete = phase(:,5) > max(phase(:,5))*(1-c_sig/100) ;
elseif abs(max(phase(:,5)))-abs(min(phase(:,5)))<-2*std(z)
rowsToDelete =  phase(:,5) < min(phase(:,5))*(1-c_sig/100);

else

    rowsToDelete = phase(:,5) > max(phase(:,5))*(1-c_sig/100)| phase(:,5) < min(phase(:,5))*(1-c_sig/100);
end
% 删除这些行
phase(rowsToDelete, :) = [];
x=phase(:,1);
xp=phase(:,2);
y=phase(:,3);
yp=phase(:,4);
z=phase(:,5);
p=phase(:,6);


j=1;


if beamratio<1

    data=[x,xp,y,yp,z,p];
total_rows = size(data, 1);
% 生成随机的行索引
num_rows_to_extract = round(total_rows*beamratio);
random_indices = randperm(total_rows, num_rows_to_extract);

% 根据随机索引从原始数据中选择行
data = data(random_indices, :);

phase=data;
x=phase(:,1);
xp=phase(:,2);
y=phase(:,3);
yp=phase(:,4);
z=phase(:,5);
p=phase(:,6);

np=size(data,1);
end

switch plotimag


    case {'exb3','exb'}



% 处理数据
x = phase(:,1);
xp = phase(:,2);
y = phase(:,3);
yp = phase(:,4);
z = phase(:,5);
p = phase(:,6);

% 转换为6D矩阵格式 [x,xp,y,yp,z,p]
X_orig = [x'; xp'; y'; yp'; z'; p'];

% 计算原始分布的Twiss参数和发射度
% X平面
sigma_x = std(X_orig(1,:));
sigma_xp = std(X_orig(2,:));
sigma_xxp = mean((X_orig(1,:)-mean(X_orig(1,:))).*(X_orig(2,:)-mean(X_orig(2,:))));
emit_x = sqrt(sigma_x^2 * sigma_xp^2 - sigma_xxp^2);
beta_x = sigma_x^2 / emit_x;
alpha_x = -sigma_xxp / emit_x;

% Y平面
sigma_y = std(X_orig(3,:));
sigma_yp = std(X_orig(4,:));
sigma_yyp = mean((X_orig(3,:)-mean(X_orig(3,:))).*(X_orig(4,:)-mean(X_orig(4,:))));
emit_y = sqrt(sigma_y^2 * sigma_yp^2 - sigma_yyp^2);
beta_y = sigma_y^2 / emit_y;
alpha_y = -sigma_yyp / emit_y;

% 输出原始Twiss参数
fprintf('\n原始分布 Twiss参数:\n');
fprintf('X平面:\n');
fprintf('发射度(x) = %.6e m·rad\n', emit_x);
fprintf('beta_x = %.6f m\n', beta_x);
fprintf('alpha_x = %.6f\n', alpha_x);
fprintf('\nY平面:\n');
fprintf('发射度(y) = %.6e m·rad\n', emit_y);
fprintf('beta_y = %.6f m\n', beta_y);
fprintf('alpha_y = %.6f\n', alpha_y);

% =========================================================================
% 核心优化区：6D 关联核密度扩束 (完美保留所有6维相关性，无需人为缩放)
% =========================================================================
N_orig = size(X_orig, 2);
N_new = N_orig * n_beam;

% 1. 计算原始 6D 分布的均值和协方差矩阵
mu_orig = mean(X_orig, 2);
Sigma_orig = cov(X_orig');

% 2. 计算 6维空间的最优平滑核宽 h (Silverman's rule of thumb)
d = 6;
h = (4 / ((d + 2) * N_orig)) ^ (1 / (d + 4)); 
tuning_factor = 0.8; % 调节因子，防止分布过度平滑
h = h * tuning_factor;

% 3. 对协方差矩阵进行 Cholesky 分解 (加极小值防奇异矩阵报错)
L = chol(Sigma_orig + eye(6)*1e-15, 'lower');

% 4. 随机抽取基准粒子 (有放回重采样)
idx = randi(N_orig, 1, N_new);
X_base = X_orig(:, idx);

% 5. 生成与原始 6D 分布具有相同相关性（相同倾角、色散）的高斯噪声
correlated_noise = L * randn(6, N_new);

% 6. 合成新粒子并进行理论严格缩放
% 缩放系数 1/sqrt(1+h^2) 保证了添加噪声后，整体的协方差(即发射度和Twiss)绝对守恒
scale_factor_6d = 1 / sqrt(1 + h^2);
X_new = mu_orig + (X_base - mu_orig + h * correlated_noise) * scale_factor_6d;

% =========================================================================

% 计算新分布的Twiss参数
sigma_x = std(X_new(1,:));
sigma_xp = std(X_new(2,:));
sigma_xxp = mean((X_new(1,:)-mean(X_new(1,:))).*(X_new(2,:)-mean(X_new(2,:))));
emit_x_new = sqrt(sigma_x^2 * sigma_xp^2 - sigma_xxp^2);
beta_x_new = sigma_x^2 / emit_x_new;
alpha_x_new = -sigma_xxp / emit_x_new;

sigma_y = std(X_new(3,:));
sigma_yp = std(X_new(4,:));
sigma_yyp = mean((X_new(3,:)-mean(X_new(3,:))).*(X_new(4,:)-mean(X_new(4,:))));
emit_y_new = sqrt(sigma_y^2 * sigma_yp^2 - sigma_yyp^2);
beta_y_new = sigma_y^2 / emit_y_new;
alpha_y_new = -sigma_yyp / emit_y_new;

% 输出新Twiss参数
fprintf('\n新分布 Twiss参数:\n');
fprintf('X平面:\n');
fprintf('发射度(x) = %.6e m·rad\n', emit_x_new);
fprintf('beta_x = %.6f m\n', beta_x_new);
fprintf('alpha_x = %.6f\n', alpha_x_new);
fprintf('\nY平面:\n');
fprintf('发射度(y) = %.6e m·rad\n', emit_y_new);
fprintf('beta_y = %.6f m\n', beta_y_new);
fprintf('alpha_y = %.6f\n', alpha_y_new);

% 保存扩展后的数据
fid = fopen(outputfile, 'w');

% 注意：修复了原代码中遇到 [] 和 nc(浮点数) 使用 %d 会报错的问题
a = {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';np*n_beam '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    if isempty(a{i,1})
        fprintf(fid, '         %s\n', a{i,2});
    else
        fprintf(fid, '%g         %s\n', a{i,1}, a{i,2}); % 使用 %g 兼容整数与浮点电荷
    end
end

X_new1 = X_new';
for j = 1:size(X_new1,1)
    fprintf(fid, '%.15f                  ', X_new1(j,:));
    fprintf(fid, '\n');
end

fprintf(fid, 'END_BUNCH \n');
fclose(fid);

% 原始分布（上排）
subplot(2,3,1);
plot(1000*X_orig(1,:), 1000*X_orig(2,:), '.b', 'MarkerSize', 1);
title('Original');
xlabel('x[mm]'); ylabel('xp[mrad]');

subplot(2,3,2);
plot(1000*X_orig(3,:), 1000*X_orig(4,:), '.b', 'MarkerSize', 1);
title('Original');
xlabel('y[mm]'); ylabel('yp[mrad]');

subplot(2,3,3);
plot(1000*X_orig(5,:), 1000*X_orig(6,:), '.b', 'MarkerSize', 1);
title('Original');
xlabel('z[mm]'); ylabel('p');

% 新分布（下排）
subplot(2,3,4);
plot(1000*X_new(1,:), 1000*X_new(2,:), '.r', 'MarkerSize', 1);
title('Expanded');
xlabel('x[mm]'); ylabel('xp[mrad]');

subplot(2,3,5);
plot(1000*X_new(3,:), 1000*X_new(4,:), '.r', 'MarkerSize', 1);
title('Expanded');
xlabel('y[mm]'); ylabel('yp[mrad]');

subplot(2,3,6);
plot(1000*X_new(5,:), 1000*X_new(6,:), '.r', 'MarkerSize', 1);
title('Expanded');
xlabel('z[mm]'); ylabel('p');



    case 'longmat'
z=z-mean(z);
        if zp0==0
[~, sorted_indices] = sort(abs(z)); % 按距离 z_mean 的绝对值排序
closest_indices = sorted_indices(1:50); % 取前 50 个粒子的索引
pz_closest_mean = mean(p(closest_indices)); % 对这 50 个粒子的 pz 求平均值
fprintf('z_mean 附近的 50 个粒子的平均 pz 值为: %.5f\n', pz_closest_mean);

% 平移 pz
p = p - pz_closest_mean;
p0=0;
        else 
p=p;
[~, sorted_indices] = sort(abs(z)); % 按距离 z_mean 的绝对值排序
closest_indices = sorted_indices(1:50); % 取前 50 个粒子的索引
pz_closest_mean = mean(p(closest_indices)); % 对这 50 个粒子的 pz 求平均值

p0=pz_closest_mean;

        end
        if fit==1
% 生成示例数据（如果您有自己的数据，请替换这部分）
x =z;
y = p;
zp=p;
% 使用多项式拟合
degree = 10; % 拟合多项式的度数
p = polyfit(x, y, degree);

% 创建拟合曲线的点
x_fit = linspace(min(x), max(x), 200);
y_fit = polyval(p, x_fit);



% 创建左侧子图

% 绘制原始数据点
plot(x, y, 'k.', 'MarkerSize', 10);
hold on;

% 绘制拟合曲线
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% 设置图形属性
xlabel('z');
ylabel('zp');
grid off;
legend('Data', 'Fitted Curve', 'Location', 'best');

% 创建标题，显示前三次项的系数
fprintf('C0=%.2f \n a1=%.2f \n a2=%.2f \n a3=%.2f \n', p(end),p(end-1),p(end-2),p(end-3));
title('initial phase space');

% 调整图形以适应标题
axis tight;
ylim_current = ylim;
ylim([ylim_current(1), ylim_current(2) + (ylim_current(2)-ylim_current(1))*0.1]);

hold off;

krf=2*pi/(3*10^8/freq/10^6);

mer1=z;
mer2=(zp+1)*energy;
 linac1=mer1;
 linac2=mer2;
linac2=mer2+EA*cos(krf.*mer1-phi_0*pi/180);
mean_ene=energy*(1+p0)+EA*cos(phi_0*pi/180);
phase(:,6)=(linac2-mean_ene)/mean_ene;
phase(:,5)=mer1;
phase(:,2)=phase(:,2)./(mean_ene).*energy;
phase(:,4)=phase(:,4)./(mean_ene).*energy;
linac22=phase(:,6);
p=phase(:,6);
z=mer1+R560*linac22+T5660*linac22.*linac22+U56660*linac22.*linac22.*linac22;
phase(:,5)=z;
mer1=z;
fid = fopen(outputfile, 'w');
a= {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';size(phase,1) '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end


for j=1:1:size(phase,1)

 
    fprintf(fid, '%.15f                 ', phase(j,:));
    fprintf(fid, '\n');
    
   
    
    
end

fprintf(fid, 'END_BUNCH \n');
fclose(fid);

z=mer1+R56*linac22+T566*linac22.*linac22+U5666*linac22.*linac22.*linac22;


zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
figure;
   
yyaxis left




[values, centers] = hist3([1000*z, phase(:,6)], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256; % 颜色映射的大小
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed]; % 在映射开始处添加白色

% 应用自定义颜色映射
colormap(blueToRed);
yyaxis left;
% 绘制色温图
imagesc(centers{:}, percentValues');  % 使用 centers 来指定 x 和 y 轴的坐标
axis xy;  % 将原点设置到左下角，并使 y 轴正向朝上
h = colorbar; % 显示颜色条

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6); % 创建6个刻度
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

xlabel('z[mm]');
ylabel('\delta');
set(gca, 'YColor', 'k');  % 确保 y 轴的颜色是黑色




    
yyaxis right
plot(1000*zhist(:,1),zhist(:,2),'k','LineWidth',2);
ylabel('current[A]');
 set(gca, 'ycolor', 'k');

        else
        
        x =z;
y = p;
zp=p;
% 使用多项式拟合



% 创建左侧子图



krf=2*pi/(3*10^8/freq/10^6);

mer1=z;
mer2=(zp+1)*energy;
 linac1=mer1;
 linac2=mer2;
linac2=mer2+EA*cos(krf.*mer1-phi_0*pi/180);
mean_ene=energy*(1+p0)+EA*cos(phi_0*pi/180);
phase(:,6)=(linac2-mean_ene)/mean_ene;
phase(:,5)=mer1;
phase(:,2)=phase(:,2)./(mean_ene).*energy;
phase(:,4)=phase(:,4)./(mean_ene).*energy;
linac22=phase(:,6);
z=mer1+R560*linac22+T5660*linac22.*linac22+U56660*linac22.*linac22.*linac22;
phase(:,5)=z;
mer1=z;
fid = fopen(outputfile, 'w');
a= {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';size(phase,1) '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end


for j=1:1:size(phase,1)

 
    fprintf(fid, '%.15f                 ', phase(j,:));
    fprintf(fid, '\n');
    
   
    
    
end

fprintf(fid, 'END_BUNCH \n');
fclose(fid);

z=phase(:,5)+R56*linac22+T566*linac22.*linac22+U5666*linac22.*linac22.*linac22;

zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];

   
yyaxis left




[values, centers] = hist3([1000*z, phase(:,6)], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256; % 颜色映射的大小
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed]; % 在映射开始处添加白色

% 应用自定义颜色映射
colormap(blueToRed);
yyaxis left;
% 绘制色温图
imagesc(centers{:}, percentValues');  % 使用 centers 来指定 x 和 y 轴的坐标
axis xy;  % 将原点设置到左下角，并使 y 轴正向朝上
h = colorbar; % 显示颜色条

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6); % 创建6个刻度
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

xlabel('z[mm]');
ylabel('\delta');
set(gca, 'YColor', 'k');  % 确保 y 轴的颜色是黑色




    
yyaxis right
plot(1000*zhist(:,1),zhist(:,2),'k','LineWidth',2);
ylabel('current[A]');
 set(gca, 'ycolor', 'k');
        
        end

sigz=std(z);

fprintf('压缩后长度为 %.5f 微米\n', sigz*10^6);

    case 'skew'

        phase=[x,xp,y,yp,z,p];


% 初始化用于存储偏度的向量
skewnessValues = zeros(1, 6);

% 计算每一维的偏度
for i = 1:6
    skewnessValues(i) = skewness(phase(:, i));
end

% 显示偏度值
fprintf('束团偏度为: ');  
fprintf('%f ', skewnessValues);

    case 'inform'

        phase=[x,xp,y,yp,z,p];

% 初始化用于存储偏度的向量
datapart = zeros(4, 6);

for i = 1:6
     datapart(1,i) = mean(phase(:, i));
end

for i = 1:6
     datapart(2,i) = std(phase(:, i));
end

% 计算每一维的偏度
for i = 1:6
    datapart(3,i) = skewness(phase(:, i));
end

for i = 1:6
     datapart(4,i) = kurtosis(phase(:, i))-3;
end

% 显示偏度值
fprintf('mean: ');  
fprintf('%f     ', datapart(1,:));
fprintf('\n');

fprintf('std: ');  
fprintf('%f     ', datapart(2,:));
fprintf('\n');

fprintf('skew: ');  
fprintf('%f     ', datapart(3,:));
fprintf('\n');

fprintf('kurtosis: ');  
fprintf('%f     ', datapart(4,:));
fprintf('\n');

phase(:,1)=phase(:,1)-mean(phase(:,1));
phase(:,2)=phase(:,2)-mean(phase(:,2));
phase(:,3)=phase(:,3)-mean(phase(:,3));
phase(:,4)=phase(:,4)-mean(phase(:,4));

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);
rms_sigxp=mean(sigxp);
rms_sigxxp=mean(sigxxp);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;

sigy=phase(:,3).*phase(:,3);
sigyp=phase(:,4).*phase(:,4);
sigyyp=phase(:,3).*phase(:,4);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig36=phase(:,3).*phase(:,6);
sig46=phase(:,4).*phase(:,6);
rms_sigy=mean(sigy);
rms_sigyp=mean(sigyp);
rms_sigyyp=mean(sigyyp);

emity=real(sqrt(rms_sigy*rms_sigyp-rms_sigyyp^2));
betyall=rms_sigy/emity;
alpyall=-rms_sigyyp/emity;

% 计算洛伦兹因子（需要定义能量）

gam=energy/0.511;
% 计算归一化发射度
emitx_norm = emitx * 1 * gam;
emity_norm = emity * 1 * gam;

fprintf('\n=== 束流参数 ===\n');
fprintf('beginning[beta_a]=%f\n',betxall);
fprintf('beginning[alpha_a]=%f\n',alpxall);
fprintf('beginning[beta_b]=%f\n', betyall);
fprintf('beginning[alpha_b]=%f\n', alpyall);

fprintf('\n=== 发射度信息 ===\n');
fprintf('X方向几何发射度: %.6e m·rad\n', emitx);
fprintf('Y方向几何发射度: %.6e m·rad\n', emity);
fprintf('X方向归一化发射度: %.6e m·rad\n', emitx_norm);
fprintf('Y方向归一化发射度: %.6e m·rad\n', emity_norm);






    case 'wgd'
e0=energy/.511;
phaseb=[x xp y yp z p];
np=size(x,1);

for j=1:1:np
z1(j)  =phaseb(j,5);
p1(j)  =phaseb(j,6)*e0+e0;
yp1(j)=phaseb(j,4);
x1(j)  =phaseb(j,1);
xp1(j)=phaseb(j,2);
y1(j)  =phaseb(j,3);
end
phase=[x1' xp1' y1' yp1' z1' p1'];
a={'?   VERSION = ' '1.0'; '?   CHARGE = ' nc; '?   COLUMNS    X    XPRIME    Y    YPRIME    Z    ' 'gamma'};

fid=fopen(outputfile, 'w');


    fprintf(fid, '%s', a{1,1});
      fprintf(fid, '%s\n', a{1,2});
      
    fprintf(fid, '%s', a{2,1});
      fprintf(fid, '%d\n', a{2,2});
      
    fprintf(fid, '%s', a{3,1});
      fprintf(fid, '%s\n', a{3,2});


fclose(fid);

fid1=fopen(outputfile, 'a');

for i=1:1:np
   fprintf(fid1, '%s   ', phase(i,:));

    fprintf(fid1, '\n');
end
% fprintf(fid1,'%s  %s  %s  %s  %s  %s \n',phase);
 fclose(fid1);


    case 'we'

        
        Data = [x xp y yp z p];

Data(:,5)=-1*(Data(:,5))/299792458;
Data(:,6)=energy/0.511*(Data(:,6)+1);





fid = fopen(outputfile, 'w');

fprintf(fid, 'SDDS1 \n');
fprintf(fid, '&parameter name=Charge, type=double, units=C, description="total charge in Coulombs" &end \n');
fprintf(fid, '&column name=x,  type=double, units=m, description="x in meters" &end \n');
fprintf(fid, '&column name=xp, type=double, description="px/pz" &end \n');
fprintf(fid, '&column name=y,  type=double, units=m, description="y in meters" &end \n');
fprintf(fid, '&column name=yp, type=double, description="py/pz" &end \n');
fprintf(fid, '&column name=t,  type=double, units=s, description="time in seconds" &end \n');
fprintf(fid, '&column name=p,  type=double, units="m$be$nc", description="relativistic gamma*beta" &end \n');
fprintf(fid, '&column name=particleID, type=ulong64, &end \n');
fprintf(fid, '&data mode=ascii &end \n');
fprintf(fid, '%s \n',nc);
fprintf(fid, '%d \n',size(z,1)-1);

for j=1:1:size(z)-1

 
    fprintf(fid, '%.30f                 ', Data(j,:));
    fprintf(fid, '%d',j);
    fprintf(fid, '\n');
    
   
    
    
end
fclose(fid);



    case 'x'

[values, centers] = hist3([1000*x, 1000*xp], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 上方的子图
subplot(10,2,[1,2]);
xmax=max(x);
xmin=min(x);
histogram(1000*x, 'BinWidth', 1000*(xmax-xmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('x');  
ylabel('Numbers');


% 下方的子图
subplot(10,2,[19,20]);
xpmax=max(xp);
xpmin=min(xp);
histogram(1000*xp, 'BinWidth', 1000*(xpmax-xpmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('xp');  
ylabel('Numbers');


% 中间的子图（粒子束密度分布）
ax2=subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
% if xrange > yrange
%     % 如果 x 范围更大，调整 y 轴限制
%     ymid = mean(ylims);
%     ylim([ymid - xrange/2, ymid + xrange/2]);
% else
%     % 如果 y 范围更大，调整 x 轴限制
%     xmid = mean(xlims);
%     xlim([xmid - yrange/2, xmid + yrange/2]);
% end
% 
% % 设置纵横比为 1:1
% pbaspect([1 1 1]);

% 其他设置
xlabel('x [mm]');
ylabel('xp [mrad]');


 case 'yz'

[values, centers] = hist3([1000*y, 1000*z], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 上方的子图
subplot(10,2,[1,2]);
xmax=max(x);
xmin=min(x);
histogram(1000*x, 'BinWidth', 1000*(xmax-xmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('x');  
ylabel('Numbers');


% 下方的子图
subplot(10,2,[19,20]);
xpmax=max(xp);
xpmin=min(xp);
histogram(1000*xp, 'BinWidth', 1000*(xpmax-xpmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('xp');  
ylabel('Numbers');


% 中间的子图（粒子束密度分布）
ax2=subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
% if xrange > yrange
%     % 如果 x 范围更大，调整 y 轴限制
%     ymid = mean(ylims);
%     ylim([ymid - xrange/2, ymid + xrange/2]);
% else
%     % 如果 y 范围更大，调整 x 轴限制
%     xmid = mean(xlims);
%     xlim([xmid - yrange/2, xmid + yrange/2]);
% end
% 
% % 设置纵横比为 1:1
% pbaspect([1 1 1]);

% 其他设置
xlabel('y [mm]');
ylabel('z [mm]');




 case 'xpyp'

[values, centers] = hist3([1000*xp, 1000*yp], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 上方的子图
subplot(10,2,[1,2]);
xmax=max(x);
xmin=min(x);
histogram(1000*x, 'BinWidth', 1000*(xmax-xmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('x');  
ylabel('Numbers');


% 下方的子图
subplot(10,2,[19,20]);
xpmax=max(xp);
xpmin=min(xp);
histogram(1000*xp, 'BinWidth', 1000*(xpmax-xpmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('xp');  
ylabel('Numbers');


% 中间的子图（粒子束密度分布）
ax2=subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
% if xrange > yrange
%     % 如果 x 范围更大，调整 y 轴限制
%     ymid = mean(ylims);
%     ylim([ymid - xrange/2, ymid + xrange/2]);
% else
%     % 如果 y 范围更大，调整 x 轴限制
%     xmid = mean(xlims);
%     xlim([xmid - yrange/2, xmid + yrange/2]);
% end
% 
% % 设置纵横比为 1:1
% pbaspect([1 1 1]);

% 其他设置
xlabel('xp [mrad]');
ylabel('yp [mrad]');






 case 'xz'



[values, centers] = hist3([1000*z, 1000*x], [tempnx tempny]);  % 修改为z和x

% 上方的子图
subplot(10,2,[1,2]);
zmax=max(z);
zmin=min(z);
histogram(1000*z, 'BinWidth', 1000*(zmax-zmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('z [mm]');  % 修改标签
ylabel('Numbers');


% 中间的子图（粒子束密度分布）
subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
if xrange > yrange
    % 如果 x 范围更大，调整 y 轴限制
    ymid = mean(ylims);
    ylim([ymid - xrange/2, ymid + xrange/2]);
else
    % 如果 y 范围更大，调整 x 轴限制
    xmid = mean(xlims);
    xlim([xmid - yrange/2, xmid + yrange/2]);
end

% 设置纵横比为 1:1
pbaspect([1 1 1]);

% 其他设置
xlabel('z [mm]');  % 修改标签
ylabel('x [mm]');  % 修改标签

% 下方的子图
subplot(10,2,[19,20]);
xmax=max(x);
xmin=min(x);
histogram(1000*x, 'BinWidth', 1000*(xmax-xmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('x [mm]');  % 修改标签
ylabel('Numbers');




    case 'xy'
        




[values, centers] = hist3([1000*x, 1000*y], [tempnx tempny]);  % 改回x和y

% 上方的子图
subplot(10,2,[1,2]);
xmax=max(x);
xmin=min(x);
histogram(1000*x, 'BinWidth', 1000*(xmax-xmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('x [mm]');  % 改回x标签
ylabel('Numbers');


% 中间的子图（粒子束密度分布）
ax2=subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
if xrange > yrange
    % 如果 x 范围更大，调整 y 轴限制
    ymid = mean(ylims);
    ylim([ymid - xrange/2, ymid + xrange/2]);
else
    % 如果 y 范围更大，调整 x 轴限制
    xmid = mean(xlims);
    xlim([xmid - yrange/2, xmid + yrange/2]);
end

% 设置纵横比为 1:1
pbaspect([1 1 1]);

% 其他设置
xlabel('x [mm]');  % 改回x标签
ylabel('y [mm]');  % 改为y标签


% 下方的子图
subplot(10,2,[19,20]);
ymax=max(y);
ymin=min(y);
histogram(1000*y, 'BinWidth', 1000*(ymax-ymin)/nbin,'FaceColor',[1,0,0]);  
xlabel('y [mm]');  % 改为y标签
ylabel('Numbers');
box on;



    case 'y'



[values, centers] = hist3([1000*y, 1000*yp], [tempnx tempny]);  % 改为y和yp

subplot(10,2,[19,20]);
ypmax=max(yp);
ypmin=min(yp);
histogram(1000*yp, 'BinWidth', 1000*(ypmax-ypmin)/nbin,'FaceColor',[1,0,0]);  
xlabel('yp [mrad]');  % 改为yp标签
ylabel('Numbers');

% 上方的子图
subplot(10,2,[1,2]);
ymax=max(y);
ymin=min(y);
histogram(1000*y, 'BinWidth', 1000*(ymax-ymin)/nbin,'FaceColor',[1,0,0]);  
xlabel('y [mm]');  % 改为y标签
ylabel('Numbers');

% 下方的子图


% 中间的子图（粒子束密度分布）
ax2=subplot(10,2,[5:16]);

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256;
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed];

% 应用自定义颜色映射
colormap(blueToRed);

% 绘制色温图
imagesc(centers{:}, percentValues');
axis xy;
h = colorbar;

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6);
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

% 获取当前的 x 和 y 轴限制
xlims = xlim();
ylims = ylim();

% 计算 x 和 y 的范围
xrange = diff(xlims);
yrange = diff(ylims);

% 确定哪个范围更大，并据此调整轴限制
if xrange > yrange
    % 如果 x 范围更大，调整 y 轴限制
    ymid = mean(ylims);
    ylim([ymid - xrange/2, ymid + xrange/2]);
else
    % 如果 y 范围更大，调整 x 轴限制
    xmid = mean(xlims);
    xlim([xmid - yrange/2, xmid + yrange/2]);
end



% 其他设置
xlabel('y [mm]');  % 改为y标签
ylabel('yp [mrad]');  % 改为yp标签




  
   
    case 'z'
        z=z-mean(z);

zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
   
yyaxis left




[values, centers] = hist3([1000*z, p], [tempnx tempny]);  % 确保第一个输入是[x,y]

% 计算每个网格内点的数量占总点数的百分比
percentValues = values / np * 100;

% 创建从蓝色到红色的渐变颜色映射，密度为0的区域为白色
n = 256; % 颜色映射的大小
blueToRed = [linspace(0, 1, n/2)', zeros(n/2, 1), linspace(1, 0, n/2)'; ...
             linspace(1, 0, n/2)', zeros(n/2, 1), zeros(n/2, 1)];
blueToRed = [1 1 1; blueToRed]; % 在映射开始处添加白色

% 应用自定义颜色映射
colormap(blueToRed);
yyaxis left;
% 绘制色温图
imagesc(centers{:}, percentValues');  % 使用 centers 来指定 x 和 y 轴的坐标
axis xy;  % 将原点设置到左下角，并使 y 轴正向朝上
h = colorbar; % 显示颜色条

% 设置颜色条刻度以显示百分比
maxPercent = max(percentValues(:));
h.Ticks = linspace(0, maxPercent, 6); % 创建6个刻度
h.TickLabels = arrayfun(@(v) sprintf('%.2f%%', v), h.Ticks, 'UniformOutput', false);

xlabel('z[mm]');
ylabel('\delta_e');
set(gca, 'YColor', 'k');  % 确保 y 轴的颜色是黑色




    
yyaxis right
plot(1000*zhist(:,1),zhist(:,2),'k','LineWidth',2);
ylabel('current[A]');
 set(gca, 'ycolor', 'k');


    
    case '1x'

zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigym;']);
phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)/sigpzm-sigpzm/sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
eval(['xxp',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
eval(['xxp',num2str(j),'(i,2)=rms_sigx;']);
eval(['xxp',num2str(j),'(i,3)=rms_sigxp;']);
eval(['norm_emit',num2str(j),'(i,2)=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));']);
 eval(['norm_emit',num2str(j),'(i,3)=rms_sigpz;']);
phase=phase1;
end


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2);%-sigxpm;
phase(:,1)=phase(:,1);%-sigxm;
phase(:,6)=phase(:,6)-sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);



emitall=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
n=size(zhist,1);
emal=ones(n,1).*emitall;
yyaxis left;
plot(1000*zhist(:,1),zhist(:,2),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('current[A]');
 set(gca, 'ycolor', 'r');
yyaxis right
plot(1000*zhist(:,1),10^6*norm_emit1(:,2),'b','LineWidth',2);
hold on;
plot(1000*zhist(:,1), 10^6*emal(:,1),'b','LineWidth',2,'LineStyle','--');
hold off;
ylabel('slice emittance[μmrad]');
 set(gca, 'ycolor', 'b');
%fprintf('sige=%f\n',rms_sigpz);

xlim([min(1000*zhist(:,1)),max(1000*zhist(:,1))]);


    case '1nx'




zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigym;']);
phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)/sigpzm-sigpzm/sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx)-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp)-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp)-mean(sig16)*mean(sig26)/mean(sigpz);
eval(['xxp',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
eval(['xxp',num2str(j),'(i,2)=rms_sigx;']);
eval(['xxp',num2str(j),'(i,3)=rms_sigxp;']);
eval(['norm_emit',num2str(j),'(i,2)=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));']);
 eval(['norm_emit',num2str(j),'(i,3)=rms_sigpz;']);
phase=phase1;
end


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)-sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx)-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp)-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp)-mean(sig16)*mean(sig26)/mean(sigpz);



emitall=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
n=size(zhist,1);
emal=ones(n,1).*emitall;
yyaxis left;
plot(1000*zhist(:,1),zhist(:,2),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('current[A]');
 set(gca, 'ycolor', 'r');
yyaxis right
plot(1000*zhist(:,1),10^6*norm_emit1(:,2),'b','LineWidth',2);
hold on;
plot(1000*zhist(:,1), 10^6*emal(:,1),'b','LineWidth',2,'LineStyle','--');
hold off;
ylabel('slice emittance[μmrad]');
 set(gca, 'ycolor', 'b');
%fprintf('sige=%f\n',rms_sigpz);

xlim([min(1000*zhist(:,1)),max(1000*zhist(:,1))]);



    case 'mex'


% 预处理
%phase(:,1) = phase(:,1) - mean(phase(:,1));
%phase(:,2) = phase(:,2) - mean(phase(:,2));
particles = [phase(:,1), phase(:,2)];

% 只保留 x > 0 的粒子
 particles = particles(particles(:,1) > 0, :);

% 按 x 坐标排序
particles = sortrows(particles, 1, 'ascend');

n = size(particles, 1);

% 预分配内存
crossSum = 0;
na=0;
% 使用向量化操作计算叉积和
for g = 1:n-1
    pi = particles(g, :);
    remainingParticles = particles(g+1:n, :);
    crossProducts = pi(1) * remainingParticles(:, 2) - pi(2) * remainingParticles(:, 1);
    na=na+size(crossProducts,1);
    crossSum = crossSum + sum(crossProducts);
end

% 计算均值发射度
emal = 2 * crossSum / n^2;

fprintf('the mean emit is %e \n', emal);


% yyaxis left;
% plot(1000*zhist(:,1),zhist(:,2),'r','LineWidth',2)
% xlabel('z[mm]');
% ylabel('current[A]');
%  set(gca, 'ycolor', 'r');
% yyaxis right
% plot(1000*zhist(:,1),meanemit(:,1),'b','LineWidth',2);
% hold on;
% plot(1000*zhist(:,1), emal(:,1),'b','LineWidth',2,'LineStyle','--');
% hold off;
% 
% 
% ylabel('slice mean emittance');
%  set(gca, 'ycolor', 'b');
% 
% 
% xlim([min(1000*zhist(:,1)),max(1000*zhist(:,1))]);



    case 'se'


        sigpzm = mean(phase(:, 6));



phaseALL(:, 6) = (phase(:, 6) - sigpzm) / (1 + sigpzm);

sigpz = phaseALL(:, 6) .* phaseALL(:, 6);
rms_sigpz = sqrt(mean(sigpz));

emitall = rms_sigpz;

fprintf('project energy spread=%f\n', emitall);

%% 去除切片中心能量偏移计算切片能散
zmin = min(z);
zmax = max(z);
binwid = (zmax - zmin) / nbin;
binEdges = linspace(zmin, zmax, nbin + 1);

% 计算每个bin的中心位置
for i = 1:nbin
    zbin(i) = zmin + (i - 0.5) * binwid;
end

% 计算电流分布
[counts, ~] = histcounts(z, binEdges);
binCounts = counts';
zhist = [zbin', binCounts ./ sum(binCounts) .* (nc) ./ (binwid / 299792458)];

% 保存原始相空间数据
phase1 = [x xp y yp z p];

% 创建新束团数据
new_x = x;
new_xp = xp;
new_y = y;
new_yp = yp;
new_z = z;
new_p = zeros(size(p));

% 按3*nbin切片处理中心能量偏移
process_nbin = 3 * nbin;
process_binwid = (zmax - zmin) / process_nbin;
process_binEdges = linspace(zmin, zmax, process_nbin + 1);

for i = 1:process_nbin
    % 找到当前切片内的粒子
    slice_mask = (z >= process_binEdges(i)) & (z < process_binEdges(i+1));
    
    if sum(slice_mask) > 0
        % 计算该切片的平均能量偏差
        slice_p = p(slice_mask);
        slice_center_p = mean(slice_p);
        
        % 减去切片中心能量
        new_p(slice_mask) = slice_p - slice_center_p;
    end
end

% 新束团相空间数据
new_phase = [new_x new_xp new_y new_yp new_z new_p];

% 现在按原始nbin计算切片能散
norm_emit1 = zeros(nbin, 3);

for i = 1:nbin
    % 找到当前切片内的粒子
    slice_mask = (new_z >= binEdges(i)) & (new_z < binEdges(i+1));
    phase_slice = new_phase(slice_mask, :);
    
    if size(phase_slice, 1) > 0
        % 计算切片中心
        sigpzm = mean(phase_slice(:, 6));
        
        % 能量转换
        emean = energy * (1 + sigpzm);
        ee = (1 + phase_slice(:, 6)) * energy;
        
        % 中心化处理
        phase_slice(:, 6) = (ee - emean) / emean;
        
        % 计算能散
        sigpz = phase_slice(:, 6) .* phase_slice(:, 6);
        rms_sigpz_slice = sqrt(mean(sigpz));
        
        % 存储结果
        norm_emit1(i, 1) = zbin(i);
        norm_emit1(i, 2) = 0;
        norm_emit1(i, 3) = rms_sigpz_slice;
    else
        norm_emit1(i, :) = [zbin(i), 0, 0];
    end
end

% 绘图
% 设置Cambridge字体
set(0, 'DefaultAxesFontName', 'Cambria');
set(0, 'DefaultTextFontName', 'Cambria');

% 绘制电流分布
yyaxis left;
plot(1000 * zhist(:, 1), zhist(:, 2), 'r', 'LineWidth', 2);
xlabel('z [mm]', 'FontName', 'Cambria', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Current [A]', 'FontName', 'Cambria', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'ycolor', 'r');
set(gca, 'FontName', 'Cambria', 'FontSize', 10);

% 绘制切片能散
yyaxis right;
plot(1000 * norm_emit1(:, 1), norm_emit1(:, 3), 'b', 'LineWidth', 2);
ylabel('Slice Energy Spread', 'FontName', 'Cambria', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'ycolor', 'b');
set(gca, 'FontName', 'Cambria', 'FontSize', 10);

% 设置图形属性
xlim([min(1000 * zhist(:, 1)), max(1000 * zhist(:, 1))]);
grid off;
set(gca, 'GridAlpha', 0.3);
box on;
set(gca, 'LineWidth', 1.5);

   case '1y'

zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[y yp x xp z p];
phase=phase1;

%sigxmall=mean(y);
%sigymall=mean(x);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigym;']);
phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)/sigpzm-sigpzm/sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
eval(['xxp',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
eval(['xxp',num2str(j),'(i,2)=rms_sigx;']);
eval(['xxp',num2str(j),'(i,3)=rms_sigxp;']);
eval(['norm_emit',num2str(j),'(i,2)=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));']);
 eval(['norm_emit',num2str(j),'(i,3)=rms_sigpz;']);
phase=phase1;
end


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)-sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);



emitall=real(energy/0.511*sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
n=size(zhist,1);
emal=ones(n,1).*emitall;
yyaxis left;
plot(1000*zhist(:,1),zhist(:,2),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('current[A]');
 set(gca, 'ycolor', 'r');
yyaxis right
plot(1000*zhist(:,1),10^6*norm_emit1(:,2),'b','LineWidth',2);
hold on;
plot(1000*zhist(:,1), 10^6*emal(:,1),'b','LineWidth',2,'LineStyle','--');
hold off;
ylabel('slice emittance[μmrad]');
 set(gca, 'ycolor', 'b');
fprintf('sige=%f\n',rms_sigpz);

xlim([min(1000*zhist(:,1)),max(1000*zhist(:,1))]);




    case '3x'




zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];

phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigym;']);
phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)/sigpzm-sigpzm/sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
eval(['xxp',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);


emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betx(i)=rms_sigx/emitx;
alpx(i)=-rms_sigxxp/emitx;

phase=phase1;
end

phase(:,1)=phase(:,1)-mean(phase(:,1));
phase(:,2)=phase(:,2)-mean(phase(:,2));
phase(:,3)=phase(:,3)-mean(phase(:,3));
phase(:,4)=phase(:,4)-mean(phase(:,4));


sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;


sigy=phase(:,3).*phase(:,3);
sigyp=phase(:,4).*phase(:,4);
sigyyp=phase(:,3).*phase(:,4);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig36=phase(:,3).*phase(:,6);
sig46=phase(:,4).*phase(:,6);
rms_sigy=mean(sigy);%-mean(sig36)*mean(sig36)/mean(sigpz);
rms_sigyp=mean(sigyp);%-mean(sig46)*mean(sig46)/mean(sigpz);
rms_sigyyp=mean(sigyyp);%-mean(sig36)*mean(sig46)/mean(sigpz);


emity=real(sqrt(rms_sigy*rms_sigyp-rms_sigyyp^2));
betyall=rms_sigy/emity;
alpyall=-rms_sigyyp/emity;





fprintf('beginning[beta_a]=%f\n',betxall);
fprintf('beginning[alpha_a]=%f\n',alpxall);
fprintf('beginning[beta_b]=%f\n', betyall);
fprintf('beginning[alpha_b]=%f\n', alpyall);


sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx)-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp)-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp)-mean(sig16)*mean(sig26)/mean(sigpz);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;

fprintf('(bmad)beginning[beta_a]=%f\n',betxall);
fprintf('(bmad)beginning[alpha_a]=%f\n',alpxall);


yyaxis left;
plot(1000*xxp1(:,1),betx(:),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('betax[m]');
 set(gca, 'ycolor', 'r')
yyaxis right
plot(1000*xxp1(:,1),alpx(:),'b','LineWidth',2);
ylabel('alpx');
 set(gca, 'ycolor', 'b')


    case'pemit'

        phase=[x xp y yp z p];
        
        phase(:,1)=phase(:,1)-mean(phase(:,1));
        phase(:,2)=phase(:,2)-mean(phase(:,2));

        phase1=phase;
        sigx=phase(:,1).*phase(:,1);
        sigxp=phase(:,2).*phase(:,2);
        sigxxp=phase(:,1).*phase(:,2);

        rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
        rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
        rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
        
        emitx0=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
        betx=rms_sigx/emitx0;
        alpx=-rms_sigxxp/emitx0;
        gamx = (1 + alpx^2) / betx;


x =phase(:,1);   % 假设的位置数据
xp = phase(:,2);  % 假设的动量数据

% 计算所有点的椭圆方程值
ellipse_values = gamx * x.^2 + 2 * alpx * x .* xp + betx * xp.^2;

% 发射度阈值
emittances = 0.25*emitx0:0.5*emitx0:maxr*emitx0;

emitmat=zeros(length(emittances),2);

% 遍历每个发射度阈值
for i = 1:length(emittances)
    % 当前发射度阈值
    emitx = emittances(i);
    
    % 找出在当前发射度椭圆内的点
    indices = ellipse_values <= emitx;
    points_inside = [x(indices), xp(indices)];
    
    % 计算这些点的实际发射度
    % 实际发射度可以通过点在相空间的分布计算得到，例如使用RMS发射度
    mean_x = mean(points_inside(:,1));
    mean_xp = mean(points_inside(:,2));
    sigma_xx = var(points_inside(:,1));
    sigma_xpxp = var(points_inside(:,2));
    sigma_xxpxp = cov(points_inside(:,1), points_inside(:,2));
    
    % 计算 RMS 发射度
    epsilon_rms = sqrt(sigma_xx * sigma_xpxp - sigma_xxpxp(1,2)^2);
    
   emitmat(i,1)=(i-1)*.5+0.25;
   emitmat(i,2)=epsilon_rms/emitx0;
end

        plot(emitmat(:,1),emitmat(:,2),'Color','r','LineWidth',2,'Marker','square');
        hold on;
        xlabel('area of ellipse/pro.emit of complete bunch');
        ylabel('pro.emit in phase ellipse/pro.emit of complete bunch');

         phase=[y yp y yp z p];
        
        phase(:,1)=phase(:,1)-mean(phase(:,1));
        phase(:,2)=phase(:,2)-mean(phase(:,2));

        phase1=phase;
        sigx=phase(:,1).*phase(:,1);
        sigxp=phase(:,2).*phase(:,2);
        sigxxp=phase(:,1).*phase(:,2);

        rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
        rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
        rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
        
        emitx0=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
        betx=rms_sigx/emitx0;
        alpx=-rms_sigxxp/emitx0;
        gamx = (1 + alpx^2) / betx;


x =phase(:,1);   % 假设的位置数据
xp = phase(:,2);  % 假设的动量数据

% 计算所有点的椭圆方程值
ellipse_values = gamx * x.^2 + 2 * alpx * x .* xp + betx * xp.^2;

% 发射度阈值
emittances = 0.25*emitx0:0.5*emitx0:maxr*emitx0;

emitmat=zeros(length(emittances),2);

% 遍历每个发射度阈值
for i = 1:length(emittances)
    % 当前发射度阈值
    emitx = emittances(i);
    
    % 找出在当前发射度椭圆内的点
    indices = ellipse_values <= emitx;
    points_inside = [x(indices), xp(indices)];
    
    % 计算这些点的实际发射度
    % 实际发射度可以通过点在相空间的分布计算得到，例如使用RMS发射度
    mean_x = mean(points_inside(:,1));
    mean_xp = mean(points_inside(:,2));
    sigma_xx = var(points_inside(:,1));
    sigma_xpxp = var(points_inside(:,2));
    sigma_xxpxp = cov(points_inside(:,1), points_inside(:,2));
    
    % 计算 RMS 发射度
    epsilon_rms = sqrt(sigma_xx * sigma_xpxp - sigma_xxpxp(1,2)^2);
    
   emitmat(i,1)=(i-1)*.5+0.25;
   emitmat(i,2)=epsilon_rms/emitx0;
end

        plot(emitmat(:,1),emitmat(:,2),'Color','b','LineWidth',2,'Marker','square');

        hold off;

        legend('xphase','yphase');

        

    case '3y'




zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];

phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigym;']);
phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;
phase(:,6)=phase(:,6)/sigpzm-sigpzm/sigpzm;
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);
eval(['xxp',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);


emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betx(i)=rms_sigx/emitx;
alpx(i)=-rms_sigxxp/emitx;

phase=phase1;
end




yyaxis left;
plot(1000*xxp1(:,1),betx(:),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('betay[m]');
 set(gca, 'ycolor', 'r')
yyaxis right
plot(1000*xxp1(:,1),alpx(:),'b','LineWidth',2);
ylabel('alpy');
 set(gca, 'ycolor', 'b')



    case '2x'

        zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigxpm;']);

phase=phase1;
end

yyaxis left;
plot(1000*zhist(:,1),center1(:,1),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('det x[m]');
 set(gca, 'ycolor', 'r')
yyaxis right
plot(1000*zhist(:,1),center1(:,2),'b','LineWidth',2);
ylabel('det xp');
 set(gca, 'ycolor', 'b')



    case '2xc'

        zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;
z0=[zhist(:,2)];
sigxmall=mean(x);
sigxpmall=mean(xp);
j=1;
for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1))-sigxmall;
sigxpm=mean(phase(:,2))-sigxpmall;


eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigxpm;']);

phase=phase1;
end

b0=[center1(:,1) center1(:,2) zhist(:,2)];
z1=zhist(:,1)-min(zhist(:,1));

for i=1:1:nbin+1

    if zhist(i,2)<0
    else
color = [z1(i)/(max(z1(:))), 0, 1 - z1(i)/(max(z1(:)))];
scatter(b0(i,1),b0(i,2),50*(z0(i)+1)/max(z0),'MarkerEdgeColor',color,'LineWidth',2);

hold on;
    end


end
box on;
xlabel('x_-x_0', 'Interpreter', 'latex');
ylabel('x''-x''_0', 'Interpreter', 'latex');


nbar=100;

Imax=max(z1);
% legend(['z=' num2str(Imax) 'A'],['z=' num2str(Imax/2) 'A'],['z=' num2str(Imax*0) 'A']);
redToBlue = [linspace(0, 1, nbar)', zeros(nbar, 1), linspace(1, 0, nbar)'];
colormap(redToBlue);
c = colorbar;
caxis([1 nbar]); % 设置颜色条的范围

c.Ticks = [1,  nbar]; % 设置颜色条的刻度
c.TickLabels = {['tail'], ['head']};


hold off;

  case '2y'

        zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;

%sigxmall=mean(x);
%sigymall=mean(y);

for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,3));
sigxpm=mean(phase(:,4));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));

eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigxpm;']);

phase=phase1;
end

yyaxis left;
plot(1000*zhist(:,1),center1(:,1),'r','LineWidth',2)
xlabel('z[mm]');
ylabel('det y[m]');
 set(gca, 'ycolor', 'r')
yyaxis right
plot(1000*zhist(:,1),center1(:,2),'b','LineWidth',2);
ylabel('det yp');
 set(gca, 'ycolor', 'b')

    case '2yc'
  zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[y yp y yp z p];
phase=phase1;
z0=[zhist(:,2)];
sigxmall=mean(y);
sigxpmall=mean(yp);
j=1;
for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];
eval(['norm_emit',num2str(j),'(i,1)=zmin+(i-1/2)*binwid;']);
n=size(phase,1);
sigxm=mean(phase(:,1))-sigxmall;
sigxpm=mean(phase(:,2))-sigxpmall;


eval(['center',num2str(j),'(i,1)=sigxm;']);
eval(['center',num2str(j),'(i,2)=sigxpm;']);

phase=phase1;
end

b0=[center1(:,1) center1(:,2) zhist(:,2)];
z1=zhist(:,1)-min(zhist(:,1));

for i=1:1:nbin+1

    if zhist(i,2)<0
    else
color = [z1(i)/(max(z1(:))), 0, 1 - z1(i)/(max(z1(:)))];
scatter(b0(i,1),b0(i,2),(z0(i)+1)*4,'MarkerEdgeColor',color,'LineWidth',2);

hold on;
    end


end
box on;

xlabel('y_{slice}-y_{mean}');
ylabel('yp_{slice}-yp_{mean}');


nbar=100;

Imax=max(z1);
% legend(['z=' num2str(Imax) 'A'],['z=' num2str(Imax/2) 'A'],['z=' num2str(Imax*0) 'A']);
redToBlue = [linspace(0, 1, nbar)', zeros(nbar, 1), linspace(1, 0, nbar)'];
colormap(redToBlue);
c = colorbar;
caxis([1 nbar]); % 设置颜色条的范围

c.Ticks = [1,  nbar]; % 设置颜色条的刻度
c.TickLabels = {['tail'], ['head']};


hold off;


    case '4x'


    zmin=min(z);
zmax=max(z);


n=size(z,1);
if round(n/1000)>200
    nbin=200;
else
    nbin=round(n/1000);
end

binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[x xp y yp z p];
phase=phase1;
Imax=max(zhist(:,2));
%Imax1=max(zhist1(:,1));
%sigxmall=mean(x);
%sigymall=mean(y);




for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];

n=size(phase,1);



sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);



emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betx=rms_sigx/emitx;
alpx=-rms_sigxxp/emitx;
gamx=(1+alpx^2)/betx;

theta = linspace(0, 2*pi, 500);
xe = sqrt(emitx * betx) * cos(theta)+sigxm;
xpe = -sqrt(emitx / betx) * (sin(theta) + alpx * cos(theta))+sigxpm;
color = [zhist(i,2)/(Imax), 0, 1 - zhist(i,2)/(Imax)];



plot(1000*xe, 1000*xpe, 'Color', color,'LineWidth',exp(2*(zhist(i,2)/Imax-.5)));
%ezplot('gamx*x^2+2*alpx*x*y+betx*y^2=emitx');
hold on;
scatter(sigxm,sigxpm,'MarkerEdgeColor',color );
hold on;
phase=phase1;
end

sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);

emitxall=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitxall;
alpxall=-rms_sigxxp/emitxall;
gamxall=(1+alpxall^2)/betxall;
theta = linspace(0, 2*pi, 500);
xeall = sqrt(emitxall * betxall) * cos(theta);
xpeall = -sqrt(emitxall / betxall) * (sin(theta) + alpxall * cos(theta));

plot(1000*xeall, 1000*xpeall, 'Color', 'k', 'LineStyle','--','LineWidth',5);

xlabel('x [mm]');
ylabel('x'' [mrad]');
%title('Phase Space Ellipse');

hold off;
nbar=100;
%legend(['I=' num2str(Imax) 'A'],['I=' num2str(Imax/2) 'A'],['I=' num2str(Imax*0) 'A']);
redToBlue = [linspace(0, 1, nbar)', zeros(nbar, 1), linspace(1, 0, nbar)'];
colormap(redToBlue);
c = colorbar;
caxis([1 nbar]); % 设置颜色条的范围

c.Ticks = [1, nbar/2, nbar]; % 设置颜色条的刻度
c.TickLabels = {[num2str(Imax*0) 'A'], [num2str(Imax/2) 'A'], [num2str(Imax) 'A']}; % 设置颜色条的标签




% sigyp=phase(:,4).*phase(:,4);
% sigyyp=phase(:,3).*phase(:,4);
% sigpz=phase(:,6).*phase(:,6);
% rms_sigpz=sqrt(mean(sigpz));
% sig36=phase(:,3).*phase(:,6);
% sig46=phase(:,4).*phase(:,6);
% rms_sigy=mean(sigy);%-mean(sig36)*mean(sig36)/mean(sigpz);
% rms_sigyp=mean(sigyp);%-mean(sig46)*mean(sig46)/mean(sigpz);
% rms_sigyyp=mean(sigyyp);%-mean(sig36)*mean(sig46)/mean(sigpz);
% 
% 
% emity=real(sqrt(rms_sigy*rms_sigyp-rms_sigyyp^2));
% betyall=rms_sigy/emity;
% alpyall=-rms_sigyyp/emity;


    case 'mma'
maxk=1;

phase=[x xp y yp z p zeros(size(z,1),1)];

zmin=min(z);
zmax=max(z);
%nbin=50;
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
[maxVal, maxRow] = max(zhist(:,2));

maxiz=zhist(maxRow,1);
minz=min(z);
maxz=max(z);

minzl=minz-maxiz;
maxzl=maxz-maxiz;

for i=1:size(phase,1)
    if phase(i,5)>maxiz
        phase(i,7)=(phase(i,5)-maxiz)/maxzl;
    else
        phase(i,7)=(phase(i,5)-maxiz)/minzl;
    end

end
phase(phase(:,7)>maxk,:)=[];

phase(:,1)=phase(:,1)-mean(phase(:,1));
phase(:,2)=phase(:,2)-mean(phase(:,2));

x=phase(:,1);
xp=phase(:,2);

z=phase(:,5);
sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);

rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;












beta_source =betxall ; % 替换为实际值
alpha_source =alpxall; % 替换为实际值
gamma_source = (1 + alpha_source^2) / beta_source;


 target_betx=1;
 target_alpx=0;


% 出口Twiss参数
beta_target = target_betx; % 替换为实际值
alpha_target = target_alpx; % 替换为实际值
gamma_target = (1 + alpha_target^2) / beta_target;
% 假设这是您的switch语句的一部分

        % 已知入口和出口的Twiss参数
        beta_s =beta_source; % 替换为实际值
        alpha_s =alpha_source; % 替换为实际值
        gamma_s =gamma_source;
        beta_t = beta_target; % 替换为实际值
        alpha_t =alpha_target; % 替换为实际值
        gamma_t = gamma_target;

        % 定义方程组
        twissEquations = @(x) [
            x(1)^2 * beta_s - 2 * x(1) * x(2) * alpha_s + x(2)^2 * gamma_s - beta_t;
            -x(1) * x(3) * beta_s + (x(1) * x(4) + x(2) * x(3)) * alpha_s - x(2) * x(4) * gamma_s - alpha_t;
            x(3)^2 * beta_s - 2 * x(3) * x(4) * alpha_s + x(4)^2 * gamma_s - gamma_t;
            x(1) * x(4) - x(2) * x(3) - 1
        ];

        % 使用fsolve求解
        M0 = [1, 0, 0, 1]; % 初始猜测
       options = optimoptions('fsolve', 'Display', 'iter','MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
M_vals = fsolve(twissEquations, M0, options);

        % 构建传输矩阵
        Mx = [M_vals(1), M_vals(2); M_vals(3), M_vals(4)];


    


xnew =[Mx(1,1), Mx(1,2)]*[x';xp'];
xpnew=[Mx(2,1), Mx(2,2)]*[x';xp'];




phase=[xnew' xpnew' phase(:,3:6)];

nc=nc*size(phase,1)/np;
x=phase(:,1);
xp=phase(:,2);
np=size(phase,1);
zmin=min(z);
zmax=max(z);
nbin=ceil(nbin*maxk);
binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;
zbin=[];
zhist=[];
for i=1:1:nbin+1
    zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];

phase1=phase;
phase=phase1;
Imax=max(zhist(:,2));
%Imax1=max(zhist1(:,1));
%sigxmall=mean(x);
%sigymall=mean(y);

rms_mis=zeros(nbin+1,1);


for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];

n=size(phase,1);



sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));



phase(:,2)=phase(:,2);
phase(:,1)=phase(:,1);

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);

rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);



emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betx=rms_sigx/emitx;
alpx=-rms_sigxxp/emitx;
gamx=(1+alpx^2)/betx;


color = [zhist(i,2)/(Imax), 0, 1 - zhist(i,2)/(Imax)];
misfact=betx*gamma_t+beta_t*gamx-2*alpha_t*alpx;
rms_mis(i)=1/2*abs(misfact);



scatter(betx-beta_t,alpx-0,50,'MarkerEdgeColor',color,'LineWidth',4 );
hold on;
phase=phase1;
end


axis equal;

xlabel('\beta_{x}-\beta_{mean}');
ylabel('\alpha_{x}-\alpha_{mean}');

rms_mis(nbin+1,:)=[];
rms_mis(1,:)=[];



rmsmis=sum(rms_mis(:,1));

r =rmsmis; % 定义圆的半径
theta = linspace(0, 2*pi, 100); % 生成角度向量

% 计算圆上的点

%title(['r = ' num2str(r)]);
hold off;
nbar=100;
%legend(['I=' num2str(Imax) 'A'],['I=' num2str(Imax/2) 'A'],['I=' num2str(Imax*0) 'A']);
redToBlue = [linspace(0, 1, nbar)', zeros(nbar, 1), linspace(1, 0, nbar)'];
colormap(redToBlue);
c = colorbar;
caxis([1 nbar]); % 设置颜色条的范围

c.Ticks = [1, nbar/2, nbar]; % 设置颜色条的刻度
c.TickLabels = {[num2str(Imax*0) 'A'], [num2str(Imax/2) 'A'], [num2str(Imax) 'A']}; % 设置颜色条的标签

ylim([-4,4]);
xlim([-4,4]);



 case 'mmz'

%% 简化的MMF计算代码

% 基本参数设置
e = 1.602e-19;
maxk = 1;

%% 1. 数据准备和分bin
phase = [x, xp, y, yp, z, p, zeros(size(z,1),1)];

zmin = min(z);
zmax = max(z);
binwid = (zmax - zmin) / nbin;
binEdges = zmin:binwid:zmax;

% 计算bin中心和电流分布
for i = 1:nbin+1
    zbin(i) = zmin + (i - 1/2) * binwid;
end

[counts, ~] = histc(z, binEdges);
zhist = [zbin', counts(:) ./ sum(counts(:)) .* nc ./ (binwid / 3 / (10^8))];

% 找到峰值电流位置
[~, maxRow] = max(zhist(:,2));
maxiz = zhist(maxRow, 1);

%% 2. 数据裁剪（可选，根据maxk参数）
minzl = zmin - maxiz;
maxzl = zmax - maxiz;

for i = 1:size(phase,1)
    if phase(i,5) > maxiz
        phase(i,7) = (phase(i,5) - maxiz) / maxzl;
    else
        phase(i,7) = (phase(i,5) - maxiz) / minzl;
    end
end
phase(phase(:,7) > maxk, :) = [];

%% 3. 计算全局X方向Twiss参数（参考1：全局）
phase(:,1) = phase(:,1) - mean(phase(:,1));
phase(:,2) = phase(:,2) - mean(phase(:,2));

rms_sigx_global = mean(phase(:,1).^2);
rms_sigxp_global = mean(phase(:,2).^2);
rms_sigxxp_global = mean(phase(:,1) .* phase(:,2));

emitx_global = real(sqrt(rms_sigx_global * rms_sigxp_global - rms_sigxxp_global^2));
betx_global = rms_sigx_global / emitx_global;
alpx_global = -rms_sigxxp_global / emitx_global;
gamx_global = (1 + alpx_global^2) / betx_global;

fprintf('全局X方向: β=%.4e m, α=%.4f, γ=%.4f, ε=%.4e m·rad\n', ...
    betx_global, alpx_global, gamx_global, emitx_global);

%% 4. 计算全局Y方向Twiss参数（参考1：全局）
phase(:,3) = phase(:,3) - mean(phase(:,3));
phase(:,4) = phase(:,4) - mean(phase(:,4));

rms_sigy_global = mean(phase(:,3).^2);
rms_sigyp_global = mean(phase(:,4).^2);
rms_sigyyp_global = mean(phase(:,3) .* phase(:,4));

emity_global = real(sqrt(rms_sigy_global * rms_sigyp_global - rms_sigyyp_global^2));
bety_global = rms_sigy_global / emity_global;
alpy_global = -rms_sigyyp_global / emity_global;
gamy_global = (1 + alpy_global^2) / bety_global;

fprintf('全局Y方向: β=%.4e m, α=%.4f, γ=%.4f, ε=%.4e m·rad\n', ...
    bety_global, alpy_global, gamy_global, emity_global);

%% 5. 找到峰值切片并计算其Twiss参数（参考2：峰值切片）
z = phase(:,5);
phase_peak = phase(z > (zbin(maxRow) - binwid/2) & z <= (zbin(maxRow) + binwid/2), :);

% X方向峰值切片
rms_sigx_peak = mean(phase_peak(:,1).^2);
rms_sigxp_peak = mean(phase_peak(:,2).^2);
rms_sigxxp_peak = mean(phase_peak(:,1) .* phase_peak(:,2));

emitx_peak = real(sqrt(rms_sigx_peak * rms_sigxp_peak - rms_sigxxp_peak^2));
betx_peak = rms_sigx_peak / emitx_peak;
alpx_peak = -rms_sigxxp_peak / emitx_peak;
gamx_peak = (1 + alpx_peak^2) / betx_peak;

fprintf('峰值切片X方向: β=%.4e m, α=%.4f, γ=%.4f, ε=%.4e m·rad\n', ...
    betx_peak, alpx_peak, gamx_peak, emitx_peak);

% Y方向峰值切片
rms_sigy_peak = mean(phase_peak(:,3).^2);
rms_sigyp_peak = mean(phase_peak(:,4).^2);
rms_sigyyp_peak = mean(phase_peak(:,3) .* phase_peak(:,4));

emity_peak = real(sqrt(rms_sigy_peak * rms_sigyp_peak - rms_sigyyp_peak^2));
bety_peak = rms_sigy_peak / emity_peak;
alpy_peak = -rms_sigyyp_peak / emity_peak;
gamy_peak = (1 + alpy_peak^2) / bety_peak;

fprintf('峰值切片Y方向: β=%.4e m, α=%.4f, γ=%.4f, ε=%.4e m·rad\n', ...
    bety_peak, alpy_peak, gamy_peak, emity_peak);

%% 6. 计算每个切片的MMF（以全局为参考）
mmf_x_global = [];
mmf_y_global = [];
zbin_center_valid = [];
current_valid = [];

phase1 = phase;

for i = 1:nbin+1
    % 选择当前切片
    phase = phase1(phase1(:,5) > (zmin + (i-1)*binwid) & ...
                   phase1(:,5) <= (zmin + i*binwid), :);
    
    if size(phase, 1) < 10  % 跳过粒子数太少的切片
        continue;
    end
    
    zbin_center_valid = [zbin_center_valid; zmin + (i - 1/2) * binwid];
    current_valid = [current_valid; zhist(i, 2)];
    
    % X方向Twiss
    rms_sigx = mean(phase(:,1).^2);
    rms_sigxp = mean(phase(:,2).^2);
    rms_sigxxp = mean(phase(:,1) .* phase(:,2));
    
    emitx = real(sqrt(rms_sigx * rms_sigxp - rms_sigxxp^2));
    betx = rms_sigx / emitx;
    alpx = -rms_sigxxp / emitx;
    gamx = (1 + alpx^2) / betx;
    
    % 计算MMF（参考：全局）
    misfactx = 1/2 * (betx*gamx_global + betx_global*gamx - 2*alpx_global*alpx);
    mmf_x_global = [mmf_x_global; misfactx + sqrt(misfactx^2 - 1)];
    
    % Y方向Twiss
    rms_sigy = mean(phase(:,3).^2);
    rms_sigyp = mean(phase(:,4).^2);
    rms_sigyyp = mean(phase(:,3) .* phase(:,4));
    
    emity = real(sqrt(rms_sigy * rms_sigyp - rms_sigyyp^2));
    bety = rms_sigy / emity;
    alpy = -rms_sigyyp / emity;
    gamy = (1 + alpy^2) / bety;
    
    % 计算MMF（参考：全局）
    misfacty = 1/2 * (bety*gamy_global + bety_global*gamy - 2*alpy_global*alpy);
    mmf_y_global = [mmf_y_global; misfacty + sqrt(misfacty^2 - 1)];
end

%% 7. 计算每个切片的MMF（以峰值切片为参考）
mmf_x_peak = [];
mmf_y_peak = [];

for i = 1:nbin+1
    % 选择当前切片
    phase = phase1(phase1(:,5) > (zmin + (i-1)*binwid) & ...
                   phase1(:,5) <= (zmin + i*binwid), :);
    
    if size(phase, 1) < 10
        continue;
    end
    
    % X方向Twiss
    rms_sigx = mean(phase(:,1).^2);
    rms_sigxp = mean(phase(:,2).^2);
    rms_sigxxp = mean(phase(:,1) .* phase(:,2));
    
    emitx = real(sqrt(rms_sigx * rms_sigxp - rms_sigxxp^2));
    betx = rms_sigx / emitx;
    alpx = -rms_sigxxp / emitx;
    gamx = (1 + alpx^2) / betx;
    
    % 计算MMF（参考：峰值切片）
    misfactx = 1/2 * (betx*gamx_peak + betx_peak*gamx - 2*alpx_peak*alpx);
    mmf_x_peak = [mmf_x_peak; misfactx + sqrt(misfactx^2 - 1)];
    
    % Y方向Twiss
    rms_sigy = mean(phase(:,3).^2);
    rms_sigyp = mean(phase(:,4).^2);
    rms_sigyyp = mean(phase(:,3) .* phase(:,4));
    
    emity = real(sqrt(rms_sigy * rms_sigyp - rms_sigyyp^2));
    bety = rms_sigy / emity;
    alpy = -rms_sigyyp / emity;
    gamy = (1 + alpy^2) / bety;
    
    % 计算MMF（参考：峰值切片）
    misfacty = 1/2 * (bety*gamy_peak + bety_peak*gamy - 2*alpy_peak*alpy);
    mmf_y_peak = [mmf_y_peak; misfacty + sqrt(misfacty^2 - 1)];
end


yyaxis left
plot(1000*zbin_center_valid, mmf_x_global, 'r-', 'LineWidth', 2.5, 'DisplayName', 'MMF_x (Global)');
hold on;
plot(1000*zbin_center_valid, mmf_y_global, 'r--', 'LineWidth', 2.5, 'DisplayName', 'MMF_y (Global)');
hold off;
xlabel('z [mm]', 'FontSize', 12);
ylabel('MMF_{slice} (Global Reference)', 'FontSize', 12);
set(gca, 'ycolor', 'r');

yyaxis right
plot(1000*zbin_center_valid, mmf_x_peak, 'b-', 'LineWidth', 2.5, 'DisplayName', 'MMF_x (Peak)');
hold on;
plot(1000*zbin_center_valid, mmf_y_peak, 'b--', 'LineWidth', 2.5, 'DisplayName', 'MMF_y (Peak)');
hold off;
ylabel('MMF_{slice} (Peak Reference)', 'FontSize', 12);
set(gca, 'ycolor', 'b');

legend('Location', 'best', 'FontSize', 10);

fprintf('\nMMF计算完成！\n');
fprintf('参考1：全局束流 - X: β=%.4e, α=%.4f | Y: β=%.4e, α=%.4f\n', ...
    betx_global, alpx_global, bety_global, alpy_global);
fprintf('参考2：峰值切片 - X: β=%.4e, α=%.4f | Y: β=%.4e, α=%.4f\n', ...
    betx_peak, alpx_peak, bety_peak, alpy_peak);


% 计算圆上的点






    case '4y'

        zmin=min(z);
zmax=max(z);


n=size(z,1);
if round(n/1000)>200
    nbin=200;
else
    nbin=round(n/1000);
end

binwid=(zmax-zmin)/nbin;
binEdges = zmin:binwid:zmax;

for i=1:1:nbin+1
     zbin(i)=zmin+(i-1/2)*(binwid);
end

[counts, ~] = histc(z, binEdges);
binCounts = counts;

zhist=[zbin',binCounts(:)./sum(binCounts(:)).*(nc)./(binwid/3/(10^8))];
phase1=[y yp x xp z p];
phase=phase1;
Imax=max(zhist(:,2));
%Imax1=max(zhist1(:,1));
%sigxmall=mean(x);
%sigymall=mean(y);




for i=1:1:nbin+1
phase(find(phase(:,5)<=(zmin+(i-1)*binwid)|phase(:,5)>(zmin+i*binwid)),:)=[];

n=size(phase,1);



sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);



emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betx=rms_sigx/emitx;
alpx=-rms_sigxxp/emitx;
gamx=(1+alpx^2)/betx;

theta = linspace(0, 2*pi, 500);
xe = sqrt(emitx * betx) * cos(theta)+sigxm;
xpe = -sqrt(emitx / betx) * (sin(theta) + alpx * cos(theta))+sigxpm;
color = [zhist(i,2)/(Imax), 0, 1 - zhist(i,2)/(Imax)];



plot(1000*xe, 10^3*xpe, 'Color', color,'LineWidth',exp(2*(zhist(i,2)/Imax-.5)));
hold on;
phase=phase1;
end


sigxm=mean(phase(:,1));
sigxpm=mean(phase(:,2));
sigpzm=mean(phase(:,6));
sigym=mean(phase(:,3));


phase(:,2)=phase(:,2)-sigxpm;
phase(:,1)=phase(:,1)-sigxm;

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);
rms_sigpz=sqrt(mean(sigpz));
sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);%-mean(sig16)*mean(sig16)/mean(sigpz);
rms_sigxp=mean(sigxp);%-mean(sig26)*mean(sig26)/mean(sigpz);
rms_sigxxp=mean(sigxxp);%-mean(sig16)*mean(sig26)/mean(sigpz);

emitxall=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitxall;
alpxall=-rms_sigxxp/emitxall;
gamxall=(1+alpxall^2)/betxall;
theta = linspace(0, 2*pi, 500);
xeall = sqrt(emitxall * betxall) * cos(theta);
xpeall = -sqrt(emitxall / betxall) * (sin(theta) + alpxall * cos(theta));

plot(1000*xeall, 10^3*xpeall, 'Color', 'k', 'LineStyle','--','LineWidth',5);

xlabel('y [mm]');
ylabel('y^{\prime} [mrad]');
%title('Phase Space Ellipse');

hold off;
nbar=100;
%legend(['I=' num2str(Imax) 'A'],['I=' num2str(Imax/2) 'A'],['I=' num2str(Imax*0) 'A']);
redToBlue = [linspace(0, 1, nbar)', zeros(nbar, 1), linspace(1, 0, nbar)'];
colormap(redToBlue);
c = colorbar;
caxis([1 nbar]); % 设置颜色条的范围

c.Ticks = [1, nbar/2, nbar]; % 设置颜色条的刻度
c.TickLabels = {[num2str(Imax*0) 'A'], [num2str(Imax/2) 'A'], [num2str(Imax) 'A']}; % 设置颜色条的标签



    case 'twm'



phase=[x xp y yp z p];

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
sigpz=phase(:,6).*phase(:,6);

sig16=phase(:,1).*phase(:,6);
sig26=phase(:,2).*phase(:,6);
rms_sigx=mean(sigx);
rms_sigxp=mean(sigxp);
rms_sigxxp=mean(sigxxp);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;

sigy=phase(:,3).*phase(:,3);
sigyp=phase(:,4).*phase(:,4);
sigyyp=phase(:,3).*phase(:,4);
sigpz=phase(:,6).*phase(:,6);
sig36=phase(:,3).*phase(:,6);
sig46=phase(:,4).*phase(:,6);
rms_sigy=mean(sigy);
rms_sigyp=mean(sigyp);
rms_sigyyp=mean(sigyyp);

emity=real(sqrt(rms_sigy*rms_sigyp-rms_sigyyp^2));
betyall=rms_sigy/emity;
alpyall=-rms_sigyyp/emity;

fprintf('betax=%f\n',betxall);
fprintf('alphax=%f\n',alpxall);
fprintf('betay=%f\n', betyall);
fprintf('alphay=%f\n', alpyall);

% =======================================================
% X 平面 Twiss 匹配 (使用解析解替代 fsolve)
% =======================================================
beta_s = betxall; 
alpha_s = alpxall; 
beta_t = target_betx; % 假设已在外部定义
alpha_t = target_alpx; % 假设已在外部定义

% 解析构建传输矩阵 Mx
Mx = [ sqrt(beta_t/beta_s),                      0;
       (alpha_s - alpha_t)/sqrt(beta_s*beta_t),  sqrt(beta_s/beta_t) ];

xnew  = Mx(1,1)*x' + Mx(1,2)*xp';
xpnew = Mx(2,1)*x' + Mx(2,2)*xp';

% =======================================================
% Y 平面 Twiss 匹配 (使用解析解替代 fsolve)
% =======================================================
beta_s_y = betyall; 
alpha_s_y = alpyall; 
beta_t_y = target_bety; % 假设已在外部定义
alpha_t_y = target_alpy; % 假设已在外部定义

% 解析构建传输矩阵 My
My = [ sqrt(beta_t_y/beta_s_y),                        0;
       (alpha_s_y - alpha_t_y)/sqrt(beta_s_y*beta_t_y), sqrt(beta_s_y/beta_t_y) ];

ynew  = My(1,1)*y' + My(1,2)*yp';
ypnew = My(2,1)*y' + My(2,2)*yp';

% =======================================================
% 更新 Phase 矩阵并重新计算目标 Twiss
% =======================================================
phase=[xnew' xpnew' ynew' ypnew' z p];

sigx=phase(:,1).*phase(:,1);
sigxp=phase(:,2).*phase(:,2);
sigxxp=phase(:,1).*phase(:,2);
rms_sigx=mean(sigx);
rms_sigxp=mean(sigxp);
rms_sigxxp=mean(sigxxp);

emitx=real(sqrt(rms_sigx*rms_sigxp-rms_sigxxp^2));
betxall=rms_sigx/emitx;
alpxall=-rms_sigxxp/emitx;

sigy=phase(:,3).*phase(:,3);
sigyp=phase(:,4).*phase(:,4);
sigyyp=phase(:,3).*phase(:,4);
rms_sigy=mean(sigy);
rms_sigyp=mean(sigyp);
rms_sigyyp=mean(sigyyp);

emity=real(sqrt(rms_sigy*rms_sigyp-rms_sigyyp^2));
betyall=rms_sigy/emity;
alpyall=-rms_sigyyp/emity;

fprintf('newbetax=%f\n',betxall);
fprintf('newalphax=%f\n',alpxall);
fprintf('newbetay=%f\n', betyall);
fprintf('newalphay=%f\n', alpyall);

% =======================================================
% 绘图部分 (完全保留)
% =======================================================
%figure; % 建议加一句 figure 防止覆盖已有图
subplot(2, 2, 1);
scatter(1000*x, 1000*xp,'red');
title('source xphase');
xlabel('x[mm]');
ylabel('xp[mrad]');
box on;

subplot(2, 2, 2);
scatter(1000*xnew', 1000*xpnew','red');
title('target xphase');
xlabel('x[mm]');
ylabel('xp[mrad]');
box on;

subplot(2, 2, 3);
scatter(1000*y, 1000*yp,'red');
title('source yphase');
xlabel('y[mm]');
ylabel('yp[mrad]');
box on;

subplot(2, 2, 4);
scatter(1000*ynew', 1000*ypnew','red');
title('target yphase');
xlabel('y[mm]');
ylabel('yp[mrad]');
box on;

% =======================================================
% 文件输出部分 (修复了 %d 无法打印 [] 和浮点数 nc 的 Bug)
% =======================================================
np = size(x,1);

fid = fopen(outputfile, 'w');
a = {'!','ASCII::3';
     '0', '   ! ix_ele'; 
     '1', '   ! n_bunch';
     np, '   ! n_particle';
     [], 'BEGIN_BUNCH';
     [], 'Electron';
     nc, '   ! bunch_charge_tot';
     -0, '        ! z_center'; 
     0, '        ! t_center'};

% 打印前 3 行字符
for i = 1:3
    fprintf(fid, '%s%s\n', a{i,1}, a{i,2});
end

% 打印第 4 到 9 行（安全处理数字和空数组）
for i = 4:9
    if isempty(a{i,1})
        fprintf(fid, '         %s\n', a{i,2});
    else
        % 使用 %g 自动适配整数 np 和极小浮点数 nc
        fprintf(fid, '%g         %s\n', a{i,1}, a{i,2}); 
    end
end

% 打印 6 列 Phase 数据
for j = 1:np
    fprintf(fid, '%.15f                  ', phase(j,:));
    fprintf(fid, '\n');
end

fprintf(fid, 'END_BUNCH\n');
fclose(fid);

 case 'exb1'

% 假设 originalParticles 是一个 Nx6 的矩阵，表示N个宏粒子在六维空间的位置
originalParticles = [x xp y yp z p]; % 您的原始宏粒子数据
numNewSamples = (n_beam-1)*size(x,1); % 新样本数量
indices = randi([1 size(originalParticles, 1)], numNewSamples, 1);
newSamples = originalParticles(indices, :) + randn(numNewSamples, 6) * perturbationScale; % 添加小的随机扰动

% 合并样本
phase = [originalParticles; newSamples];

subplot(2,3,1);
scatter(x,xp,'red');
box on;
subplot(2,3,2);
scatter(y,yp,'red');
box on;
subplot(2,3,3);
scatter(z,p,'red');
box on;
subplot(2,3,4);
scatter(phase(:,1),phase(:,2),'red');
box on;
subplot(2,3,5);
scatter(phase(:,3),phase(:,4),'red');
box on;
subplot(2,3,6);
scatter(phase(:,5),phase(:,6),'red');
box on;

fid = fopen(outputfile, 'w');
a= {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';np*n_beam '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end


for j=1:1:size(phase,1)

 
    fprintf(fid, '%.15f                 ', phase(j,:));
    fprintf(fid, '\n');
    
   
    
    
end

fprintf(fid, 'END_BUNCH \n');
fclose(fid);

    case 'exb2'

        % 计算均值和协方差矩阵
originalPoints = [x xp y yp z p]; % 您的原始宏粒子数据

mu = mean(originalPoints);
Sigma = cov(originalPoints);

% 自举重采样
numNewSamples = (n_beam-1)*size(x,1); % 新样本数量等于原始样本数量
indices = randi([1 size(originalPoints, 1)], numNewSamples, 1);
bootstrapSamples = originalPoints(indices, :);

% 添加随机扰动
% 这里的扰动系数（perturbationFactor）可以根据需要调整

randomPerturbation = mvnrnd(zeros(1, 6), Sigma, numNewSamples) * perturbationScale;
newSamples = bootstrapSamples + randomPerturbation;
% 合并样本
phase = [originalPoints; newSamples];

subplot(2,3,1);
scatter(x,xp,'red');
title('source xphase');
box on;
subplot(2,3,2);
scatter(y,yp,'red');
title('source yphase');
box on;
subplot(2,3,3);

scatter(z,p,'red');
title('source zphase');
box on;
subplot(2,3,4);

scatter(phase(:,1),phase(:,2),'red');
title('expand xphase');
box on;
subplot(2,3,5);

scatter(phase(:,3),phase(:,4),'red');
title('expand yphase');
box on;
subplot(2,3,6);

scatter(phase(:,5),phase(:,6),'red');
title('expand zphase');
box on;

fid = fopen(outputfile, 'w');
a= {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';np*n_beam '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end


for j=1:1:size(phase,1)

 
    fprintf(fid, '%.15f                 ', phase(j,:));
    fprintf(fid, '\n');
    
   
    
    
end

fprintf(fid, 'END_BUNCH \n');
fclose(fid);


    case 'wgb'

p=energy/0.511.*p+energy/0.511;
phase=[x xp y yp z p];
numSlices = nbin;
results = []; % 改为动态数组

% 计算 z 的范围和每个切片的边界
zMin = min(phase(:,5));
zMax = max(phase(:,5));
zEdges = linspace(zMin, zMax, numSlices + 1);

% 计算切片长度
slice_length = (zMax - zMin) / numSlices;

% 循环处理每个切片
for i = 1:numSlices
    % 找到当前切片的粒子
    sliceMask = phase(:,5) >= zEdges(i) & phase(:,5) < zEdges(i+1);
    slice = phase(sliceMask, :);

    % 检查切片粒子数是否小于5
    numMacroParticles = size(slice, 1);
    if numMacroParticles < 5
        continue; % 粒子数小于5，直接跳过
    end

    % 计算切片的物理参数
    s = zMin+(i-1/2)*(zMax-zMin)/numSlices; % 切片的平均 z 位置
    gamma = mean(slice(:,6)); % 切片的平均 gamma
    dgamma = std(slice(:,6)); % gamma 的标准差

    % 计算 x 和 y 方向的归一化发射度（正确公式）
    % 中心化坐标
    x_centered = slice(:,1);
    xp_centered = slice(:,2);
    y_centered = slice(:,3);
    yp_centered = slice(:,4);
    


    x_centered1 = slice(:,1)-mean(slice(:,1));
    xp_centered1 = slice(:,2)-mean(slice(:,2));
    y_centered1 = slice(:,3)-mean(slice(:,3));
    yp_centered1 = slice(:,4)-mean(slice(:,4));
    
    % 计算 RMS 量
    x2 = mean(x_centered1.^2);
    xp2 = mean(xp_centered1.^2);
    xxp = mean(x_centered1.*xp_centered1);
    
    y2 = mean(y_centered1.^2);
    yp2 = mean(yp_centered1.^2);
    yyp = mean(y_centered1.*yp_centered1);
    
    % 几何发射度
    emit_x_geom = sqrt(x2 * xp2 - xxp^2);
    emit_y_geom = sqrt(y2 * yp2 - yyp^2);
    
    % 归一化发射度
    xemit = gamma * emit_x_geom;
    yemit = gamma * emit_y_geom;

    % 计算 x 和 y 方向的 RMS 值
    xrms = sqrt(x2);
    yrms = sqrt(y2);

    % 计算束团中心
    xavg = mean(slice(:,1));
    yavg = mean(slice(:,3));

    xpavg = gamma *mean(slice(:,2));
    ypavg = gamma *mean(slice(:,4));

    % 计算 Twiss 参数 alpha
    alphax = -xxp / emit_x_geom;
    alphay = -yyp / emit_y_geom;

    % 计算切片电流
    current = numMacroParticles / np * nc / (slice_length / 299792458);

   % fprintf('c=%f emit= %f %f \n', current, 1e7*xemit, 1e7*yemit);

    % 将计算结果添加到结果矩阵
    results = [results; s, gamma, dgamma, xemit, yemit, xrms, yrms, xavg, yavg, xpavg, ypavg, alphax, alphay, current, 0, numMacroParticles];
end

% 文件写入部分
fid = fopen(outputfile, 'w');
actualBins = size(results, 1); % 实际输出的切片数

fprintf(fid, '%d\n', actualBins); % 输出实际的切片数

for j = 1:actualBins
    fprintf(fid, '%d    ', results(j, 1:end-1));
    fprintf(fid, '\n');
end
    case 'ewb'

        nc=charge;

% 读取文本文件数据

data=[x,xp,y,yp,z,p];
total_rows = size(data, 1);
% 生成随机的行索引
num_rows_to_extract = round(total_rows*beamratio);
random_indices = randperm(total_rows, num_rows_to_extract);

% 根据随机索引从原始数据中选择行
data = data(random_indices, :);

np1=size(data,1);

np=size(data,1);

fid = fopen(outputfile, 'w');

a= {'!','ASCII::3';'0' '   ! ix_ele'; '1' '   ! n_bunch';np '   ! n_particle';[] 'BEGIN_BUNCH';[] 'Electron';nc '   ! bunch_charge_tot';-0 '        ! z_center'; 0 '        ! t_center'};
for i=1:1:3
    fprintf(fid, '%s', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end
for i=4:1:9
    fprintf(fid, '%d         ', a{i,1});
    fprintf(fid, '%s\n', a{i,2});
end


for j=1:1:np1

 
    fprintf(fid, '%.15f                 ', data(j,1:end));
    fprintf(fid, '\n');
    
   
    
    
end





fprintf(fid, 'END_BUNCH\n');

% 关闭文件
fclose(fid);






end









figHandles = findall(groot, 'Type', 'figure');

% 检查是否有图形窗口打开
if ~isempty(figHandles)

% 获取屏幕大小
% 使用你指定的弹窗大小

% 获取屏幕尺寸
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);

% 计算居中位置
center_x = (screen_width - pnx) / 2;
center_y = (screen_height - pny) / 2;

% 获取当前图形句柄并设置为指定大小和居中位置
fig = gcf;
set(fig, 'Units', 'pixels', 'Position', [center_x, center_y, pnx, pny]);
set(fig, 'Color', 'white');

% 获取所有坐标轴
ax = findall(fig, 'type', 'axes');

% 检查是否有绘图内容
hasPlotContent = false;
for i = 1:length(ax)
    % 检查每个坐标轴是否有子对象（线条、点等）
    children = get(ax(i), 'Children');
    if ~isempty(children)
        hasPlotContent = true;
        break;
    end
end

if hasPlotContent
    % 遍历所有坐标轴对象
    for i = 1:length(ax)
        % 设置当前坐标轴
        axes(ax(i));
        
        % 设置标题为 Cambria Math 粗体并调整大小（如果存在）
        title_obj = get(gca, 'Title');
        if ~isempty(title_obj.String)
            set(title_obj, 'FontName', 'Cambria Math', 'FontWeight', 'bold', 'FontSize', 14);
        end
        
        % 设置 x 轴标签为 Cambria Math 粗体并调整大小
        set(get(gca, 'XLabel'), 'FontName', 'Cambria Math', 'FontWeight', 'bold');
        
        % 检查是否使用了 yyaxis
        if strcmp(get(gca, 'YAxisLocation'), 'left')
            % 单 y 轴情况
            set(get(gca, 'YLabel'), 'FontName', 'Cambria Math', 'FontWeight', 'bold');
        else
            % 双 y 轴情况
            % 设置左侧 y 轴标签为 Cambria Math 粗体并调整大小
            yyaxis left
            set(get(gca, 'YLabel'), 'FontName', 'Cambria Math', 'FontWeight', 'bold');

            % 设置右侧 y 轴标签为 Cambria Math 粗体并调整大小
            yyaxis right
            set(get(gca, 'YLabel'), 'FontName', 'Cambria Math', 'FontWeight', 'bold');
        end
        
        % 设置图例为 Cambria Math 粗体并调整大小
        leg = findobj(ax(i), 'Type', 'Legend');
        if ~isempty(leg)
            set(leg, 'FontName', 'Cambria Math', 'FontWeight', 'bold');
        end
        
        % 设置坐标轴刻度标签字体
        set(gca, 'FontName', 'Cambria Math', 'FontWeight', 'bold');
        
        % 设置整体字体大小（适应小窗口）
        set(gca, 'FontSize', pny/35);
        
        % 设置线宽
        set(gca, 'LineWidth', 1);
        
        % 设置坐标轴颜色
        %set(gca, 'XColor', 'k', 'YColor', 'k');
        
        % 开启方框
        box on;
    end

    % 处理总标题字体（如果存在）
    sgtitle_obj = findall(fig, 'Type', 'text', 'Tag', 'sgtitle');
    if ~isempty(sgtitle_obj)
        set(sgtitle_obj, 'FontName', 'Cambria Math', 'FontWeight', 'bold', 'FontSize', 12);
    end

    % 刷新图形
    drawnow;
    
else
    % 如果没有绘图内容，在图形窗口中央显示 "Done!"
    clf; % 清除当前图形
    text(0.5, 0.5, 'Done!', 'FontName', 'Cambria Math', 'FontWeight', 'bold', ...
         'FontSize', 43, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'Units', 'normalized');
    axis off; % 隐藏坐标轴
end

else

fprintf('Done!');

end






