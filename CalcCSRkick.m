function [bend_integrals, segment_integrals, min_R16] = CalcCSRkick(beamline1, h_value,sigz ,calculate_min_R16)
    % CalcCSRkick_Fixed
    % 修复了上一版中 R_cumulative 未在弯铁后更新的 Bug，同时保留了计算加速。

    beamline=rev_beamline(beamline1);
    
    if nargin < 4, calculate_min_R16 = 0; end
    if nargin < 3, fprintf('WE NEED THE BUNCH LENGTH!!'); end
    
    % --- 预处理 ---
    [bend_rhos, bend_indices] = extractBendInfo(beamline);
    num_bends = length(bend_rhos);
    
    % 内存预分配
    min_R16 = 0;
    R16_storage = [];
    R16_counter = 0;
    if calculate_min_R16 == 1
        estimated_steps = length(beamline) * 501; 
        R16_storage = zeros(estimated_steps, 1);
    end
    
    intB_R51 = 0; intB_R52 = 0;
    sum_segment_R51 = 0; sum_segment_R52 = 0;
    
    R_cumulative = eye(6);
    
    % 段间状态
    in_segment = false;
    segment_start_R = [];
    segment_length = 0;
    current_bend_count = 0; 
    
    % --- 主循环 ---
    for elem_idx = 1:length(beamline)
        element = beamline(elem_idx);
        
        % 判断是否为弯铁
        is_bend = false;
        if current_bend_count < num_bends && elem_idx == bend_indices(current_bend_count + 1)
            is_bend = true;
        end
        
        if is_bend
            current_bend_count = current_bend_count + 1;
            
            % 1. 先处理之前的段间 (如果有)
            if in_segment && segment_length > 0
                rho_for_segment = bend_rhos(current_bend_count);
                [dR51, dR52] = calculateSegmentContribution(segment_start_R, segment_length, h_value, sigz, rho_for_segment);
                sum_segment_R51 = sum_segment_R51 + dR51;
                sum_segment_R52 = sum_segment_R52 + dR52;
            end
            
            % 2. 计算弯铁内部积分 (使用快速矩阵缓存)
            [dR51, dR52, R16_vals_bend] = integrateBendElement_Fast(element, R_cumulative, h_value, calculate_min_R16);
            intB_R51 = intB_R51 + dR51;
            intB_R52 = intB_R52 + dR52;
            
            % 存储 R16
            if calculate_min_R16 == 1 && ~isempty(R16_vals_bend)
                n_new = length(R16_vals_bend);
                if R16_counter + n_new > length(R16_storage)
                    R16_storage = [R16_storage; zeros(length(R16_storage), 1)];
                end
                R16_storage(R16_counter+1 : R16_counter+n_new) = R16_vals_bend;
                R16_counter = R16_counter + n_new;
            end
            
            % 3. === 关键修复：更新 R_cumulative ===
            % 必须使用元件的完整矩阵更新位置，否则下一个元件的起点是错的
            R_cumulative = getElementMatrix(element) * R_cumulative;
            
            % 重置段间状态
            in_segment = false;
            segment_length = 0;
            
        else
            % 非弯铁元件 (Drift, Quad 等)
            
            % 1. R16 采样与更新
            if calculate_min_R16 == 1 && element.length > 0
                % 既计算了采样，也返回了更新后的 R_cumulative
                [R_cumulative_new, R16_vals_drift] = propagateAndSampleDrift(element, R_cumulative);
                
                % 存储 R16
                n_new = length(R16_vals_drift);
                if R16_counter + n_new > length(R16_storage)
                    R16_storage = [R16_storage; zeros(length(R16_storage), 1)];
                end
                R16_storage(R16_counter+1 : R16_counter+n_new) = R16_vals_drift;
                R16_counter = R16_counter + n_new;
                
                % 更新主矩阵
                R_cumulative = R_cumulative_new;
            else
                % 不需要采样，直接传输
                R_cumulative = getElementMatrix(element) * R_cumulative;
            end
            
            % 2. 维护段间逻辑
            % 只有当后面还有弯铁时才累积段间
            if current_bend_count < num_bends
                if ~in_segment
                    % 开始新的段间
                    % 注意：segment_start_R 应该是段间的"起点"。
                    % 在上面的代码中，R_cumulative 刚刚被更新为当前元件的"出口"状态。
                    % 如果这是段里的第一个元件，段的起点应该是这个元件的"入口"状态吗？
                    % 不，原代码逻辑是：R_cumulative 随着 loop 推进。
                    % 当 in_segment 变为 true 时，记录 segment_start_R = R_cumulative。
                    % 此时 R_cumulative 是上一个元件(非弯铁)的出口。
                    % 所以这逻辑是对的。
                    
                    % 实际上，原代码是在 loop 最后更新 R。
                    % 如果是弯铁 -> 结束段间 -> 计算弯铁 -> 更新 R。
                    % 如果是非弯铁 -> 更新 R -> (如果是新段间, start=R)。
                    % 等等，原代码中：
                    % R_element = ...
                    % R_cumulative = R_element * R_cumulative (在 loop 最末尾)
                    % 这里的顺序是：先处理段间/弯铁积分，然后更新 R。
                    % 所以 segment_start_R 记录的是 *当前元件进入前* 的状态吗？
                    
                    % 让我们看原代码：
                    % if ~in_segment -> in_segment=true; segment_start_R = R_cumulative;
                    % R_cumulative = R_element * R_cumulative;
                    % 对！segment_start_R 记录的是段间 **第一个元件入口处** 的矩阵。
                    
                    % 在我的 Fixed 代码中，上面已经执行了 R_cumulative = ... * R_cumulative。
                    % 所以现在的 R_cumulative 是 **出口** 状态。
                    % 这会导致错误！我们需要回退或者记录入口状态。
                
                    % === 逻辑修正 ===
                    % 我们需要记录元件更新前的 R。
                end 
            end
        end
    end
    
    % --- 逻辑修正：重写主循环以严格匹配原代码顺序 ---
    % 上面的分叉逻辑太容易出错了，我们改回原代码的结构，只把耗时的函数替换掉。
    
    % 调用下面的 CalcCSRkick_Fixed_Strict
    [bend_integrals, segment_integrals, min_R16] = CalcCSRkick_Fixed_Strict(beamline, h_value, calculate_min_R16, sigz, bend_rhos, bend_indices);
end

function [bend_integrals, segment_integrals, min_R16] = CalcCSRkick_Fixed_Strict(beamline, h_value, calculate_min_R16, sigz, bend_rhos, bend_indices)
    % 严格遵循原版逻辑顺序的优化版
    
    min_R16 = 0;
    R16_storage = [];
    R16_counter = 0;
    if calculate_min_R16 == 1, R16_storage = zeros(length(beamline)*500, 1); end

    intB_R51 = 0; intB_R52 = 0;
    sum_segment_R51 = 0; sum_segment_R52 = 0;
    
    R_cumulative = eye(6);
    in_segment = false;
    segment_start_R = [];
    segment_length = 0;
    current_bend_index = 0;
    
    num_bends = length(bend_rhos);
    
    for elem_idx = 1:length(beamline)
        element = beamline(elem_idx);
        
        % 判断是否为弯铁 (利用索引加速)
        is_bend = (current_bend_index < num_bends) && (elem_idx == bend_indices(current_bend_index + 1));
        
        if is_bend
            current_bend_index = current_bend_index + 1;
            
            % 处理段间
            if in_segment && segment_length > 0
                rho_for_segment = bend_rhos(current_bend_index);
                [dR51, dR52] = calculateSegmentContribution(segment_start_R, segment_length, h_value, sigz, rho_for_segment);
                sum_segment_R51 = sum_segment_R51 + dR51;
                sum_segment_R52 = sum_segment_R52 + dR52;
            end
            
            % 处理弯铁积分 (使用优化后的函数)
            [dR51, dR52, R16_vals] = integrateBendElement_Fast(element, R_cumulative, h_value, calculate_min_R16);
            intB_R51 = intB_R51 + dR51;
            intB_R52 = intB_R52 + dR52;
            
            % 收集 R16
            if calculate_min_R16 == 1 && ~isempty(R16_vals)
                n = length(R16_vals);
                if R16_counter + n > length(R16_storage), R16_storage = [R16_storage; zeros(n,1)]; end
                R16_storage(R16_counter+1:R16_counter+n) = R16_vals;
                R16_counter = R16_counter + n;
            end
            
            in_segment = false;
            segment_length = 0;
            
        else
            % 非弯铁
            
            % 收集 R16 (使用优化后的采样)
            if calculate_min_R16 == 1 && element.length > 0
                % 注意：这里只采样，不更新 R_cumulative，因为 R_cumulative 统一在最后更新
                [~, R16_vals] = propagateAndSampleDrift(element, R_cumulative); 
                
                n = length(R16_vals);
                if R16_counter + n > length(R16_storage), R16_storage = [R16_storage; zeros(n,1)]; end
                R16_storage(R16_counter+1:R16_counter+n) = R16_vals;
                R16_counter = R16_counter + n;
            end
            
            if current_bend_index < num_bends
                if ~in_segment
                    in_segment = true;
                    segment_start_R = R_cumulative; % 这里的 R_cumulative 是进入当前元件前的状态，完全正确
                    segment_length = element.length;
                else
                    segment_length = segment_length + element.length;
                end
            end
        end
        
        % 统一更新传输矩阵 (恢复原版逻辑)
        % 这样保证了 drift 和 bend 的 update 逻辑与原版一致
        R_cumulative = getElementMatrix(element) * R_cumulative; 
    end
    
    bend_integrals = [intB_R51, intB_R52];
    segment_integrals = [sum_segment_R51, sum_segment_R52];
    if calculate_min_R16 == 1, min_R16 = min(R16_storage(1:R16_counter)); end
end

% --- 辅助函数 ---

function [bend_rhos, bend_indices] = extractBendInfo(beamline)
    bend_rhos = []; bend_indices = [];
    for i = 1:length(beamline)
        el = beamline(i);
        is_bend = false;
        % 逻辑同前
        if isfield(el, 'h') && abs(el.h) > 1e-10
            is_bend = true; rho = 1/el.h;
        elseif isfield(el, 'angle') && el.angle ~= 0
            if isfield(el, 'length') && el.length~=0, rho=el.length/el.angle; else, rho=1.0; end
            is_bend = true;
        elseif isfield(el, 'type')
             if any(strcmpi(el.type, {'bend', 'dipole', 'sbend', 'rbend', 'sector', 'rectangular'}))
                 is_bend = true; if isfield(el, 'h') && el.h~=0, rho = 1/el.h; else, rho = 1.0; end
             end
        end
        if is_bend
            bend_rhos = [bend_rhos; rho]; bend_indices = [bend_indices; i];
        end
    end
end

function [dR51, dR52, R16_vals] = integrateBendElement_Fast(element, R_initial, h_value, calculate_R16)
    if isfield(element, 'h'), h_bend = element.h; else, h_bend = 0; end
    n_steps = 500;
    ds = element.length / n_steps;
    
    % 缓存矩阵 (假设弯铁是均匀的)
    temp_el = element; 
    temp_el.length = ds;
    M_step = getElementMatrix(temp_el);
    
    temp_el.length = ds/2;
    M_half = getElementMatrix(temp_el);
    
    dR51 = 0; dR52 = 0;
    R16_vals = [];
    if calculate_R16 == 1, R16_vals = zeros(n_steps, 1); end
    
    if abs(h_bend) > 1e-10, csr_factor = abs(h_bend)^(2/3); else, csr_factor = 1; end
    
    R_current = R_initial;
    for i = 1:n_steps
        R_mid = M_half * R_current; % 进半步
        
        R51_mid = R_mid(5, 1); R52_mid = R_mid(5, 2); R56_mid = R_mid(5, 6);
        C_mid = 1 / (1 - h_value * R56_mid);
        val = C_mid^(4/3) * csr_factor * ds;
        
        dR51 = dR51 + R51_mid * val;
        dR52 = dR52 + R52_mid * val;
        
        if calculate_R16 == 1, R16_vals(i) = R_mid(1, 6); end
        
        R_current = M_step * R_current; % 进整步 (为下一次循环准备)
    end
end

function [dR51, dR52] = calculateSegmentContribution(segment_start_R, segment_length, h_value, sigz, rho_typical)
    R51 = segment_start_R(5, 1); R52 = segment_start_R(5, 2); R56 = segment_start_R(5, 6);
    C = 1 / (1 - h_value * R56);
    phi0 = (6 * sigz / abs(rho_typical))^(1/3);
    phi_i = phi0 / (C^(1/3));
    term = 2 * segment_length / (abs(rho_typical) * phi_i) + 1;
    if term <= 0, term = 1e-10; end
    dR51 = R51 * C * log(term);
    dR52 = R52 * C * log(term);
end

function [R_final, R16_vals] = propagateAndSampleDrift(element, R_start)
    n_samples = 10;
    R16_vals = zeros(n_samples, 1);
    ds = element.length / n_samples;
    temp_el = element; temp_el.length = ds;
    M_step = getElementMatrix(temp_el); % 缓存
    
    R_curr = R_start;
    for i = 1:n_samples
        R_curr = M_step * R_curr;
        R16_vals(i) = R_curr(1, 6);
    end
    R_final = R_curr;
end