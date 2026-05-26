clc;

% 直接指定输入文件路径
input_file = 'gamma.lte';  % 请修改为你的文件路径

% 检查文件是否存在
if ~exist(input_file, 'file')
    fprintf('错误：文件不存在: %s\n', input_file);
    return;
end

% 读取文件
file_content = fileread(input_file);

% 预处理文件内容
processed_content = preprocess_elegant_file(file_content);

% 解析元件和束线
[elements, beamlines] = parse_elegant_content(processed_content);

% 在命令行显示Bmad格式
display_bmad_lattice(elements, beamlines);

function processed = preprocess_elegant_file(content)
    % 移除注释并处理多行定义
    
    % 按行分割
    lines = strsplit(content, '\n');
    clean_lines = {};
    
    for i = 1:length(lines)
        line = strtrim(lines{i});
        
        % 跳过空行和注释行
        if isempty(line) || startsWith(line, '%') || startsWith(line, '!')
            continue;
        end
        
        % 移除行内注释
        comment_pos = strfind(line, '%');
        if ~isempty(comment_pos)
            line = strtrim(line(1:comment_pos(1)-1));
        end
        
        if ~isempty(line)
            clean_lines{end+1} = line;
        end
    end
    
    % 合并多行定义（以&结尾或逗号结尾但下一行有缩进）
    merged_lines = {};
    current_def = '';
    
    for i = 1:length(clean_lines)
        line = clean_lines{i};
        
        % 如果当前行是多行定义的延续
        if ~isempty(current_def)
            current_def = [current_def ' ' strtrim(line)];
        else
            current_def = line;
        end
        
        % 检查是否需要继续下一行
        if endsWith(strtrim(current_def), '&') || ...
           (endsWith(strtrim(current_def), ',') && i < length(clean_lines))
            % 移除行尾的&
            current_def = strrep(current_def, '&', '');
            continue;
        else
            % 定义完成
            merged_lines{end+1} = current_def;
            current_def = '';
        end
    end
    
    if ~isempty(current_def)
        merged_lines{end+1} = current_def;
    end
    
    processed = strjoin(merged_lines, '\n');
end

function [elements, beamlines] = parse_elegant_content(content)
    % 使用容器数组而不是结构体
    elements = containers.Map();
    beamlines = containers.Map();
    beamline_order = {};
    
    lines = strsplit(content, '\n');
    
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end
        
        % 查找冒号分隔
        colon_pos = strfind(line, ':');
        if isempty(colon_pos)
            continue;
        end
        
        element_name = strtrim(line(1:colon_pos(1)-1));
        definition = strtrim(line(colon_pos(1)+1:end));
        
        % 移除末尾分号
        if endsWith(definition, ';')
            definition = definition(1:end-1);
        end
        
        % 判断是束线还是元件
        if contains(upper(definition), 'LINE=')
            % 束线定义
            beamline = parse_beamline_def(definition);
            beamlines(element_name) = beamline;
            beamline_order{end+1} = element_name;
        else
            % 元件定义
            element = parse_element_def(definition);
            if ~isempty(element)
                elements(element_name) = element;
            end
        end
    end
    
    % 将beamline_order存储到beamlines中
    beamlines('_ORDER_') = beamline_order;
end

function element = parse_element_def(definition)
    element = struct();
    element.params = struct();
    
    % 按逗号分割，但要小心处理引号内的逗号
    parts = smart_split(definition, ',');
    
    % 第一部分是类型
    element.type = upper(strtrim(parts{1}));
    
    % 解析参数
    for i = 2:length(parts)
        param = strtrim(parts{i});
        if isempty(param)
            continue;
        end
        
        eq_pos = strfind(param, '=');
        if ~isempty(eq_pos)
            param_name = upper(strtrim(param(1:eq_pos(1)-1)));
            param_value = strtrim(param(eq_pos(1)+1:end));
            
            % 移除引号
            if startsWith(param_value, '"') && endsWith(param_value, '"')
                param_value = param_value(2:end-1);
            end
            
            % 解析Elegant表达式
            num_value = parse_elegant_expression(param_value);
            if ~isnan(num_value)
                element.params.(param_name) = num_value;
            else
                element.params.(param_name) = param_value;
            end
        end
    end
end

function parts = smart_split(str, delimiter)
    % 智能分割，处理引号内的分隔符
    parts = {};
    current_part = '';
    in_quotes = false;
    
    for i = 1:length(str)
        char = str(i);
        if char == '"'
            in_quotes = ~in_quotes;
            current_part = [current_part char];
        elseif char == delimiter && ~in_quotes
            parts{end+1} = current_part;
            current_part = '';
        else
            current_part = [current_part char];
        end
    end
    
    if ~isempty(current_part)
        parts{end+1} = current_part;
    end
end

function value = parse_elegant_expression(expr)
    % 解析Elegant表达式，按照逆波兰记法的规则
    % 例如："22.5 pi * 3 * 180 /" = ((22.5 * pi) * 3) / 180
    % 例如："pi 0 *" = pi * 0
    % 例如："2 pi * 180 /" = (2 * pi) / 180
    
    % 移除首尾空格
    expr = strtrim(expr);
    
    % 如果是空字符串，返回0
    if isempty(expr)
        value = 0;
        return;
    end
    
    % 如果是纯数字，直接返回
    num_value = str2double(expr);
    if ~isnan(num_value)
        value = num_value;
        return;
    end
    
    try
        % 按空格分割表达式
        tokens = strsplit(expr);
        tokens = tokens(~cellfun(@isempty, tokens)); % 移除空元素
        
        % 使用栈来计算表达式
        stack = [];
        
        for i = 1:length(tokens)
            token = tokens{i};
            
            if strcmp(token, '*')
                % 乘法：取栈顶两个数相乘
                if length(stack) >= 2
                    b = stack(end);
                    a = stack(end-1);
                    stack = stack(1:end-2);
                    stack(end+1) = a * b;
                end
            elseif strcmp(token, '/')
                % 除法：取栈顶两个数相除
                if length(stack) >= 2
                    b = stack(end);
                    a = stack(end-1);
                    stack = stack(1:end-2);
                    if b ~= 0
                        stack(end+1) = a / b;
                    else
                        stack(end+1) = 0;
                    end
                end
            elseif strcmp(token, '+')
                % 加法
                if length(stack) >= 2
                    b = stack(end);
                    a = stack(end-1);
                    stack = stack(1:end-2);
                    stack(end+1) = a + b;
                end
            elseif strcmp(token, '-')
                % 减法
                if length(stack) >= 2
                    b = stack(end);
                    a = stack(end-1);
                    stack = stack(1:end-2);
                    stack(end+1) = a - b;
                end
            else
                % 数值或常数
                if strcmp(token, 'pi')
                    stack(end+1) = pi;
                else
                    num = str2double(token);
                    if ~isnan(num)
                        stack(end+1) = num;
                    end
                end
            end
        end
        
        % 返回栈顶结果
        if ~isempty(stack)
            value = stack(end);
        else
            value = 0;
        end
        
    catch ME
        fprintf('Warning: Could not parse expression "%s": %s\n', expr, ME.message);
        value = 0;
    end
end

function beamline = parse_beamline_def(definition)
    % 提取LINE=后面的内容
    line_pos = strfind(upper(definition), 'LINE=');
    if isempty(line_pos)
        beamline = {};
        return;
    end
    
    content = strtrim(definition(line_pos+5:end));
    
    % 移除括号
    content = strrep(content, '(', '');
    content = strrep(content, ')', '');
    
    % 按逗号分割
    elements = strsplit(content, ',');
    beamline = {};
    
    for i = 1:length(elements)
        elem = strtrim(elements{i});
        if ~isempty(elem)
            beamline{end+1} = elem;
        end
    end
end

function display_bmad_lattice(elements, beamlines)
    % 在开头添加Q: marker
    fprintf('Q: marker\n\n');
    
    % 转换并显示元件
    element_names = keys(elements);
    
    for i = 1:length(element_names)
        name = element_names{i};
        element = elements(name);
        
        bmad_line = convert_to_bmad(name, element);
        if ~isempty(bmad_line)
            fprintf('%s\n', bmad_line);
        end
    end
    
    % 显示束线
    if isKey(beamlines, '_ORDER_')
        beamline_order = beamlines('_ORDER_');
        if ~isempty(beamline_order)
            fprintf('\n');
            
            for i = 1:length(beamline_order)
                name = beamline_order{i};
                beamline = beamlines(name);
                
                % 构建束线字符串
                line_str = sprintf('%s: line = (', name);
                current_line_length = length(line_str);
                
                for j = 1:length(beamline)
                    elem = beamline{j};
                    if j > 1
                        elem_str = [', ' elem];
                    else
                        elem_str = [' ' elem];
                    end
                    
                    % 检查是否需要换行（假设每行最大80字符）
                    if current_line_length + length(elem_str) > 80 && j > 1
                        fprintf('%s,&\n', line_str);
                        line_str = elem;
                        current_line_length = length(elem);
                    else
                        line_str = [line_str elem_str];
                        current_line_length = current_line_length + length(elem_str);
                    end
                end
                
                fprintf('%s)\n', line_str);
            end
            
            % use语句放在最后
            fprintf('\nuse, %s\n', beamline_order{end});
        end
    end
end

function bmad_line = convert_to_bmad(name, element)
    bmad_line = '';
    
    switch element.type
        case {'DRIFT', 'CSRDRIFT'}
            L = get_param_value(element.params, 'L', 0);
            bmad_line = sprintf('%s: drift, l = %g', name, L);
            
        case 'QUAD'
            L = get_param_value(element.params, 'L', 0);
            K1 = get_param_value(element.params, 'K1', 0);
            bmad_line = sprintf('%s: quadrupole, l = %g, k1 = %g', name, L, K1);
            
        case 'SEXT'
            L = get_param_value(element.params, 'L', 0);
            K2 = get_param_value(element.params, 'K2', 0);
            bmad_line = sprintf('%s: sextupole, l = %g, k2 = %g', name, L, K2);
            
        case {'SBEND', 'RBEND', 'CSRCSBEND'}
            L = get_param_value(element.params, 'L', 0);
            ANGLE = get_param_value(element.params, 'ANGLE', 0);
            K1 = get_param_value(element.params, 'K1', 0);
            E1 = get_param_value(element.params, 'E1', 0);
            E2 = get_param_value(element.params, 'E2', 0);
            
            bmad_line = sprintf('%s: sbend, l = %g, angle = %g', name, L, ANGLE);
            if K1 ~= 0
                bmad_line = sprintf('%s, k1 = %g', bmad_line, K1);
            end
            if E1 ~= 0
                bmad_line = sprintf('%s, e1 = %g', bmad_line, E1);
            end
            if E2 ~= 0
                bmad_line = sprintf('%s, e2 = %g', bmad_line, E2);
            end
            
        case {'CSRCSBEN', 'CSBEND', 'CSRBEND'}
            L = get_param_value(element.params, 'L', 0);
            ANGLE = get_param_value(element.params, 'ANGLE', 0);
            K1 = get_param_value(element.params, 'K1', 0);
            E1 = get_param_value(element.params, 'E1', 0);
            E2 = get_param_value(element.params, 'E2', 0);
            
            bmad_line = sprintf('%s: sbend, l = %g, angle = %g', name, L, ANGLE);
            if K1 ~= 0
                bmad_line = sprintf('%s, k1 = %g', bmad_line, K1);
            end
            if E1 ~= 0
                bmad_line = sprintf('%s, e1 = %g', bmad_line, E1);
            end
            if E2 ~= 0
                bmad_line = sprintf('%s, e2 = %g', bmad_line, E2);
            end
            
        case {'HKICKER', 'VKICKER', 'MONITOR'}
            L = get_param_value(element.params, 'L', 0);
            bmad_line = sprintf('%s: drift, l = %g', name, L);
            
        case 'MARK'
            bmad_line = sprintf('%s: marker', name);
            
        otherwise
            bmad_line = sprintf('! %s: %s (unsupported)', name, element.type);
    end
end

function value = get_param_value(params, param_name, default_value)
    if isfield(params, param_name)
        value = params.(param_name);
        if isnan(value) || isempty(value)
            value = default_value;
        end
    else
        value = default_value;
    end
end