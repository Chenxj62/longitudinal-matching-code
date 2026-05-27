function [zc, mean_delta, I_slice, mean_delta_raw] = bmad2slice( ...
    infile, nslices, fit_order, min_particles,source_file)
% BMAD2SLICE  Read a BMAD-format beam file, bin into longitudinal slices,
%             fit slice-mean delta with a high-order weighted polynomial,
%             evaluate the fitted curve at slice centers, optionally write
%             to disk, and return the slice arrays.
%
%   [z, delta_fit, I, delta_raw] = bmad2slice(infile, outfile, nslices, ...
%                                             fit_order, fit_mode, min_particles)
%
%     infile        : BMAD ASCII beam file
%     outfile       : output slice data file, 3 columns: z, delta_fit, I.
%                     Pass '' to skip file writing.
%     nslices       : number of longitudinal slices, default 100.
%     fit_order     : polynomial fitting order, default 5.
%     fit_mode      : 'poly' or 'none'. default 'poly'.
%     min_particles : drop edge slices whose particle count N is below
%                     this threshold. 0 disables. default 0.
%
%   Output columns:
%     z [m]     : slice center, longitudinal position
%     delta     : fitted slice-mean delta if fit_mode='poly';
%                 otherwise raw slice-mean delta.
%     I [A]     : slice current = q_macro * N_slice * c / dz
%
%   The polynomial fit is weighted by slice particle number N.

    if nargin < 2 || isempty(nslices)
        nslices = 100;
    end

    if nargin < 3 || isempty(fit_order)
        fit_order = 5;
    end


        fit_mode = 'poly';
   
    if nargin < 4 || isempty(min_particles)
        min_particles = 0;
    end

    % ---------- Read BMAD beam file ----------
    if source_file==1
    [np] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 3);
    [nc] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 6);

    [x, xp, y, yp, z, p] = textread(infile, ...
        '%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f', ...
        np, 'headerlines', 9); %#ok<NASGU,ASGLU>

    % ---------- Constants ----------
    c = 2.99792458e8;     % speed of light [m/s]

    % ---------- Reference momentum and delta ----------
    % If the 6th column in your BMAD file is already delta, keep this:
    delta = p;

    % If the 6th column is momentum p instead of delta, use this instead:
    % p0 = mean(p);
    % delta = (p - p0) ./ p0;

    % ---------- Macro-particle charge ----------
    Q_total     = nc;
    q_per_macro = Q_total / np;

    else
nc=source_file; %%%%%bunch charge [C]

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




% 计算第六列与平均值的相对误差
relativeError =(extractedData(:, 6) - energy) ./ energy;

energy=energy/1e6;

% 更新第六列元素
extractedData(:, 6) = 1*relativeError;

extractedData(1,:)=[];
extractedData(:,1)=extractedData(:,1)-mean(extractedData(:,1));
extractedData(:,3)=extractedData(:,3)-mean(extractedData(:,3));
extractedData(:,2)=extractedData(:,2)-mean(extractedData(:,2));
extractedData(:,4)=extractedData(:,4)-mean(extractedData(:,4));

z  = extractedData(:, 5);
p  = extractedData(:, 6);

 % ---------- Constants ----------
    c = 2.99792458e8;     % speed of light [m/s]

    % ---------- Reference momentum and delta ----------
    % If the 6th column in your BMAD file is already delta, keep this:
    delta = p;

    % If the 6th column is momentum p instead of delta, use this instead:
    % p0 = mean(p);
    % delta = (p - p0) ./ p0;

    % ---------- Macro-particle charge ----------
    Q_total     = nc;
    q_per_macro = Q_total / np;



    end

    % ---------- Slice in z ----------
    zmin  = min(z);
    zmax  = max(z);
    edges = linspace(zmin, zmax, nslices + 1);
    dz    = edges(2) - edges(1);

    zc = 0.5 * (edges(1:end-1) + edges(2:end));
    zc = zc(:);

    % Bin assignment for each particle
    bin = min(floor((z - zmin) / dz) + 1, nslices);
    bin(bin < 1) = 1;

    N         = accumarray(bin, 1,     [nslices 1]);
    sum_delta = accumarray(bin, delta, [nslices 1]);

    mean_delta_raw = zeros(nslices, 1);
    mask = N > 0;
    mean_delta_raw(mask) = sum_delta(mask) ./ N(mask);

    % ---------- Slice current ----------
    I_slice = q_per_macro .* N .* c ./ dz;

    % ---------- Trim edge bins with too few particles ----------
    n_dropped_head = 0;
    n_dropped_tail = 0;

    if min_particles > 0
        keep = N >= min_particles;
        i_first = find(keep, 1, 'first');
        i_last  = find(keep, 1, 'last');

        if isempty(i_first)
            error('bmad2slice: no slice has N >= %d particles.', min_particles);
        end

        n_dropped_head = i_first - 1;
        n_dropped_tail = nslices - i_last;

        idx_keep = i_first:i_last;

        zc             = zc(idx_keep);
        mean_delta_raw = mean_delta_raw(idx_keep);
        I_slice        = I_slice(idx_keep);
        N              = N(idx_keep);
        nslices        = numel(idx_keep);
    end

    % ---------- High-order fitted delta(z) ----------
    mean_delta = mean_delta_raw;

    switch lower(fit_mode)
        case 'none'
            mean_delta = mean_delta_raw;

        case 'poly'
            y = mean_delta_raw(:);
            w = N(:);

            valid = w > 0 & isfinite(y) & isfinite(zc);

            if nnz(valid) < fit_order + 1
                error('bmad2slice: not enough valid slices for polynomial order %d.', fit_order);
            end

            % Avoid too-high polynomial order.
            fit_order_eff = min(fit_order, nnz(valid) - 1);

            % Weighted high-order polynomial fit.
            mean_delta = weighted_polyfit_eval(zc, y, w, fit_order_eff);

        otherwise
            error('fit_mode must be ''poly'' or ''none''.');
    end

    % ---------- Write output ----------
   

    % ---------- Return as column vectors ----------
    zc             = zc(:);
    mean_delta     = mean_delta(:);
    I_slice        = I_slice(:);
    mean_delta_raw = mean_delta_raw(:);
end


function yfit = weighted_polyfit_eval(x, y, w, order)
%WEIGHTED_POLYFIT_EVAL Weighted polynomial fit and evaluation.
%
%   Fits y(x) with a polynomial of given order using weighted least squares.
%   The fitted curve is evaluated at the original x points.
%
%   Important:
%     x is internally normalized to avoid numerical ill-conditioning.

    x = x(:);
    y = y(:);
    w = w(:);

    valid = w > 0 & isfinite(x) & isfinite(y) & isfinite(w);

    xv = x(valid);
    yv = y(valid);
    wv = w(valid);

    if numel(xv) < order + 1
        error('weighted_polyfit_eval: not enough valid points.');
    end

    % Normalize x for numerical stability.
    x_mean = mean(xv);
    x_std  = std(xv);

    if x_std == 0
        error('weighted_polyfit_eval: x has zero standard deviation.');
    end

    xn  = (xv - x_mean) ./ x_std;
    xna = (x  - x_mean) ./ x_std;

    % Vandermonde matrix: [x^order, ..., x, 1]
    A = zeros(numel(xn), order + 1);
    for j = 0:order
        A(:, j+1) = xn.^(order - j);
    end

    % Weighted least squares:
    % min sum_i w_i * (A_i c - y_i)^2
    sw = sqrt(wv);
    Aw = A .* sw;
    yw = yv .* sw;

    coeff = Aw \ yw;

    % Evaluate fitted polynomial at all original x points.
    Aa = zeros(numel(xna), order + 1);
    for j = 0:order
        Aa(:, j+1) = xna.^(order - j);
    end

    yfit = Aa * coeff;
end