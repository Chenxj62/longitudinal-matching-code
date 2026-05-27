function [z, delta, Q_total, np] = bmad2mp(infile, source_file)
% BMAD2MP  Read a BMAD- or ASTRA-format beam file and return per-particle
%          longitudinal data (no slicing, no fitting).
%
%   [z, delta, Q_total, np] = bmad2mp(infile, source_file)
%
%     infile      : particle file
%     source_file : 1            -> BMAD ASCII beam file
%                   numeric (~=1) -> ASTRA file; the value is interpreted
%                                    as the bunch charge in [C]
%
%   Outputs:
%     z       [N x 1]  longitudinal position [m] (mean-centered)
%     delta   [N x 1]  relative energy deviation (p - <p>) / <p>
%     Q_total scalar   total bunch charge [C]
%     np      scalar   number of valid macro-particles
%
%   Companion to bmad2slice.m -- same I/O conventions, but returns the
%   raw particle arrays instead of slice-averaged quantities. Use this
%   to feed csr_calc_mp / run_csr_mp / lps_sim_csr_mp.

    % ---------- Read beam file ----------
    if source_file == 1
        [np] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 3);
        [nc] = textread(infile, '%f%*s%*s%*s', 1, 'headerlines', 6);

        [x, xp, y, yp, z, p] = textread(infile, ...
            '%n%n%n%n%n%n%*f%*f%*f%*f%*f%*f%*f', ...
            np, 'headerlines', 9); %#ok<NASGU,ASGLU>

        % If the 6th column in your BMAD file is already delta, keep this:
        delta = p;

        % If the 6th column is momentum p instead of delta, use this instead:
        % p0 = mean(p);
        % delta = (p - p0) ./ p0;

        Q_total = nc;

    else
        nc = source_file;     %%%%%bunch charge [C]

        % 读取文本文件数据
        data = dlmread(infile);

        % 根据条件筛选行
        filteredData = data(data(:, 10) > 0, :);

        % 提取前六列数据
        extractedData = filteredData(:, [1, 4, 2, 5, 3, 6]);

        np = size(extractedData, 1) - 1;

        % 第一列与第一行的差值
        extractedData(:, 1) = extractedData(:, 1) + extractedData(1, 1);

        % 第二列与第一行的差值
        extractedData(:, 2) = (extractedData(:, 2) + extractedData(1, 2)) / extractedData(1, 6);

        % 第四列与第一行的差值并除以第一行的第六列
        extractedData(:, 3) = extractedData(:, 3) + extractedData(1, 3);

        % 第五列与第一行的差值并除以第一行的第六列
        extractedData(:, 4) = (extractedData(:, 4) + extractedData(1, 4)) / extractedData(1, 6);

        % 计算第六列的平均值
        extractedData(:, 6) = extractedData(:, 6) + extractedData(1, 6);
        energy = mean(extractedData(:, 6));

        % 计算第六列与平均值的相对误差
        relativeError = (extractedData(:, 6) - energy) ./ energy;

        % 更新第六列元素
        extractedData(:, 6) = 1 * relativeError;

        extractedData(1, :) = [];
        extractedData(:, 1) = extractedData(:, 1) - mean(extractedData(:, 1));
        extractedData(:, 3) = extractedData(:, 3) - mean(extractedData(:, 3));
        extractedData(:, 2) = extractedData(:, 2) - mean(extractedData(:, 2));
        extractedData(:, 4) = extractedData(:, 4) - mean(extractedData(:, 4));

        z     = extractedData(:, 5);
        delta = extractedData(:, 6);

        Q_total = nc;
    end

    % ---------- Center z, return column vectors ----------
    z     = z(:);
    z     = z - mean(z);
    delta = delta(:);

    np = numel(z);
end
