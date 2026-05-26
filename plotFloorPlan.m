function result = plotFloorPlan(beamline, doPlot, options)
% plotFloorPlan - 绘制束线平面布局图 (包含一阶束流轨道叠加与放大)
%
% 输入参数:
%   beamline - 束线
%   doPlot   - 1 为画图并返回完整数据；其他值只返回边界和坐标数据
%   options  - 绘图与轨道配置结构体
%
% =========================================
% OPTIONS 功能一览表:
% =========================================
% [基础布局配置]
% options.figure_size       - 图形窗口大小，默认 [600, 400]
% options.line_width        - 参考轨道的线宽，默认 1
% options.axis_line_width   - 坐标轴的线宽，默认 1
% options.font_name         - 图窗字体，默认 'Cambria Math'
% options.font_size         - 字体大小，默认自适应窗口比例
% options.title_text        - 图像标题，默认 'Beamline Floor Plan'
% options.title_font_size   - 标题字体大小，默认自适应窗口比例
% options.rect_width_factor - 直线元件绘图宽度的乘积因子，默认 1
% options.rect_length       - 元件在图纸上的物理厚度(外径)，默认 1.7 (m)
% options.focus_elements    - 需要渲染形状的元件类型列表，默认 {'quadrupole', 'dipole', 'sextupole', 'cavity', 'TGU'}
% options.ang               - 束线整体在物理空间中的起始旋转角度(度)，默认 0
%
% [束流轨道配置]
% options.plot_orbit        - 是否推演并绘制一阶束流轨道，默认 false
% options.initial_orbit     - 初始束流状态 [x; px; y; py; z; dp]，默认 zeros(6,1)
% options.initial_gamma     - 初始相对论因子 gamma，默认 100
% options.orbit_color       - 轨道线条颜色，默认 'r' (红色)
% options.orbit_scale       - 轨道横向偏移(x)在图纸上的视觉放大倍数，默认 50
%
% [自动颜色]
% 非drift/marker/edge的未知元件类型自动检测并分配颜色渲染为矩形，
% 同一类型颜色一致，不同类型从内置调色板依次取色，无需手动指定。
%
% 输出:
% result.element_positions  - 包含每个元件的 name, type, pos_start(入口),
%                             pos_center(中点), pos_end(出口) 坐标
% ====================================================

if nargin < 2 || isempty(doPlot), doPlot = 1; end
if nargin < 3, options = struct(); end
if ~isfield(options, 'figure_size'), options.figure_size = [600, 400]; end
if ~isfield(options, 'line_width'), options.line_width = 1; end
if ~isfield(options, 'axis_line_width'), options.axis_line_width = 1; end
if ~isfield(options, 'font_name'), options.font_name = 'Cambria Math'; end
if ~isfield(options, 'font_size'), options.font_size = min(options.figure_size)/40; end
if ~isfield(options, 'title_text'), options.title_text = 'Beamline Floor Plan'; end
if ~isfield(options, 'title_font_size'), options.title_font_size = min(options.figure_size)/40; end
if ~isfield(options, 'rect_width_factor'), options.rect_width_factor = 1; end
if ~isfield(options, 'rect_length'), options.rect_length = 1.7; end
if ~isfield(options, 'focus_elements'), options.focus_elements = {'quadrupole', 'dipole', 'sextupole', 'cavity', 'TGU'}; end
if ~isfield(options, 'ang'), options.ang = 0; end
if ~isfield(options, 'plot_orbit'), options.plot_orbit = false; end
if ~isfield(options, 'initial_orbit'), options.initial_orbit = zeros(6, 1); end
if ~isfield(options, 'initial_gamma'), options.initial_gamma = 100; end
if ~isfield(options, 'orbit_color'), options.orbit_color = 'r'; end
if ~isfield(options, 'orbit_scale'), options.orbit_scale = 50; end

% 自动颜色调色板（用于未知元件类型）
auto_color_palette = [
    0.85, 0.33, 0.10;   % 橙
    0.00, 0.60, 0.60;   % 青
    0.80, 0.60, 0.00;   % 金
    0.40, 0.70, 0.30;   % 浅绿
    0.60, 0.40, 0.80;   % 淡紫
    0.90, 0.50, 0.50;   % 粉
    0.30, 0.30, 0.30;   % 深灰
    0.00, 0.45, 0.70;   % 钢蓝
    0.70, 0.20, 0.50;   % 酒红
    0.50, 0.80, 0.80;   % 浅青
];
auto_color_map = containers.Map();  % type -> color 映射
auto_color_idx = 0;

rotation_rad = options.ang * pi / 180;
max_points = length(beamline) * 25;
floorplan = struct();
floorplan.x = zeros(max_points, 1);
floorplan.y = zeros(max_points, 1);
floorplan.s = zeros(max_points, 1);
floorplan.theta = zeros(max_points, 1);
floorplan.elements = struct('name', {}, 'type', {}, 'pos_start', [], 'pos_end', [], 'x_center', [], 'y_center', [], 's_center', [], 'theta_center', []);
floorplan.element_positions = struct('name', {}, 'type', {}, 'pos_start', [], 'pos_center', [], 'pos_end', []);

current_x = 0; current_y = 0; current_s = 0; current_theta = 0;
point_index = 1;
floorplan.x(point_index) = current_x;
floorplan.y(point_index) = current_y;
floorplan.s(point_index) = current_s;
floorplan.theta(point_index) = current_theta;
point_index = point_index + 1;
element_index = 1;
all_elem_idx = 1;

if options.plot_orbit
    current_V = options.initial_orbit(:);
    current_gamma = options.initial_gamma;
    orbit_x_list = zeros(max_points, 1);
    orbit_y_list = zeros(max_points, 1);
    orbit_s_list = zeros(max_points, 1);
    orbit_state_list = zeros(max_points, 6);
    orbit_x_list(1) = current_x - (current_V(1) * options.orbit_scale) * sin(current_theta);
    orbit_y_list(1) = current_y + (current_V(1) * options.orbit_scale) * cos(current_theta);
    orbit_s_list(1) = current_s;
    orbit_state_list(1, :) = current_V';
    orbit_idx = 2;
end

for i = 1:length(beamline)
    element = beamline(i);
    if ~isfield(element, 'type'), continue; end
    if ~isfield(element, 'length'), element.length = 0; end
    if ~isfield(element, 'name'), element.name = sprintf('Unnamed_%s_%d', element.type, i); end
    if ~isfield(element, 'h'), element.h = 0; end

    element_type = element.type;
    element_length = element.length;
    element_start_x = current_x;
    element_start_y = current_y;
    element_start_theta = current_theta;
    old_point_index = point_index;

    switch element_type
        case 'dipole'
            angle = 0;
            if element.h ~= 0
                rho = 1 / element.h;
                angle = -element_length / rho;
            end
            if abs(angle) < 1e-10 || element_length <= 0
                new_x = current_x + element_length * cos(current_theta);
                new_y = current_y + element_length * sin(current_theta);
                floorplan.x(point_index) = new_x; floorplan.y(point_index) = new_y;
                floorplan.s(point_index) = current_s + element_length;
                floorplan.theta(point_index) = current_theta;
                point_index = point_index + 1;
            else
                radius = element_length / angle;
                circle_center_x = current_x - radius * sin(current_theta);
                circle_center_y = current_y + radius * cos(current_theta);
                start_angle = current_theta - pi/2;
                end_angle = start_angle + angle;
                num_samples = 20;
                angles = linspace(start_angle, end_angle, num_samples);
                for j = 1:num_samples
                    floorplan.x(point_index) = circle_center_x + radius * cos(angles(j));
                    floorplan.y(point_index) = circle_center_y + radius * sin(angles(j));
                    floorplan.s(point_index) = current_s + (j-1) * element_length / (num_samples-1);
                    floorplan.theta(point_index) = angles(j) + pi/2;
                    point_index = point_index + 1;
                end
                current_theta = current_theta + angle;
            end
        otherwise
            new_x = current_x + element_length * cos(current_theta);
            new_y = current_y + element_length * sin(current_theta);
            floorplan.x(point_index) = new_x; floorplan.y(point_index) = new_y;
            floorplan.s(point_index) = current_s + element_length;
            floorplan.theta(point_index) = current_theta;
            point_index = point_index + 1;
    end

    current_x = floorplan.x(point_index-1);
    current_y = floorplan.y(point_index-1);
    current_s = current_s + element_length;

    if options.plot_orbit
        num_points_added = point_index - old_point_index;
        if num_points_added > 1
            slice_element = element;
            slice_element.length = element_length / (num_points_added - 1);
            for j = 1:num_points_added
                if j == 1, continue; end
                R_slice = getElementMatrix(slice_element, current_gamma);
                current_V = R_slice * current_V;
                ref_x = floorplan.x(old_point_index + j - 1);
                ref_y = floorplan.y(old_point_index + j - 1);
                ref_theta = floorplan.theta(old_point_index + j - 1);
                ref_s = floorplan.s(old_point_index + j - 1);
                orbit_x_list(orbit_idx) = ref_x - (current_V(1) * options.orbit_scale) * sin(ref_theta);
                orbit_y_list(orbit_idx) = ref_y + (current_V(1) * options.orbit_scale) * cos(ref_theta);
                orbit_s_list(orbit_idx) = ref_s;
                orbit_state_list(orbit_idx, :) = current_V';
                orbit_idx = orbit_idx + 1;
            end
        else
            R = getElementMatrix(element, current_gamma);
            current_V = R * current_V;
            if strcmp(element_type, 'cavity')
                current_gamma = current_gamma / R(6,6);
            end
            ref_x = floorplan.x(point_index - 1);
            ref_y = floorplan.y(point_index - 1);
            ref_theta = floorplan.theta(point_index - 1);
            ref_s = floorplan.s(point_index - 1);
            orbit_x_list(orbit_idx) = ref_x - (current_V(1) * options.orbit_scale) * sin(ref_theta);
            orbit_y_list(orbit_idx) = ref_y + (current_V(1) * options.orbit_scale) * cos(ref_theta);
            orbit_s_list(orbit_idx) = ref_s;
            orbit_state_list(orbit_idx, :) = current_V';
            orbit_idx = orbit_idx + 1;
        end
    end

    if strcmp(element_type, 'dipole') && element.h ~= 0
        angle = -element_length * element.h;
        radius = element_length / angle;
        circle_center_x = element_start_x - radius * sin(element_start_theta);
        circle_center_y = element_start_y + radius * cos(element_start_theta);
        mid_angle = element_start_theta - pi/2 + angle/2;
        element_center_x = circle_center_x + radius * cos(mid_angle);
        element_center_y = circle_center_y + radius * sin(mid_angle);
        element_center_theta = mid_angle + pi/2;
    else
        element_center_x = (element_start_x + current_x) / 2;
        element_center_y = (element_start_y + current_y) / 2;
        element_center_theta = element_start_theta;
    end

    floorplan.element_positions(all_elem_idx).name = element.name;
    floorplan.element_positions(all_elem_idx).type = element_type;
    floorplan.element_positions(all_elem_idx).pos_start = [element_start_x, element_start_y];
    floorplan.element_positions(all_elem_idx).pos_center = [element_center_x, element_center_y];
    floorplan.element_positions(all_elem_idx).pos_end = [current_x, current_y];
    all_elem_idx = all_elem_idx + 1;

    % ====== 改动：非drift、非marker、非edge且有长度的元件全部收集 ======
    skip_types = {'drift', 'marker', 'edge'};
    should_collect = ismember(element_type, options.focus_elements) || ...
                     (~ismember(lower(element_type), skip_types) && element_length > 1e-10);

    if should_collect
        floorplan.elements(element_index).name = element.name;
        floorplan.elements(element_index).type = element_type;
        floorplan.elements(element_index).pos_start = [element_start_x, element_start_y];
        floorplan.elements(element_index).pos_end = [current_x, current_y];
        floorplan.elements(element_index).x_center = element_center_x;
        floorplan.elements(element_index).y_center = element_center_y;
        floorplan.elements(element_index).theta_center = element_center_theta;
        floorplan.elements(element_index).length = element_length;
        floorplan.elements(element_index).h = element.h;
        if strcmp(element_type, 'cavity')
            if isfield(element, 'frequency')
                floorplan.elements(element_index).num_cells = round(element.length/(1.5e8/element.frequency/1e6));
            else
                floorplan.elements(element_index).num_cells = 1;
            end
        end
        element_index = element_index + 1;
    end
end

floorplan.x = floorplan.x(1:point_index-1);
floorplan.y = floorplan.y(1:point_index-1);
floorplan.s = floorplan.s(1:point_index-1);
floorplan.theta = floorplan.theta(1:point_index-1);

if options.plot_orbit
    floorplan.orbit_x_global = orbit_x_list(1:orbit_idx-1);
    floorplan.orbit_y_global = orbit_y_list(1:orbit_idx-1);
    floorplan.orbit_s = orbit_s_list(1:orbit_idx-1);
    floorplan.orbit_state = orbit_state_list(1:orbit_idx-1, :);
end

if abs(rotation_rad) > 1e-10
    x_rot = floorplan.x * cos(rotation_rad) - floorplan.y * sin(rotation_rad);
    y_rot = floorplan.x * sin(rotation_rad) + floorplan.y * cos(rotation_rad);
    floorplan.x = x_rot; floorplan.y = y_rot;

    if options.plot_orbit
        ox_rot = floorplan.orbit_x_global * cos(rotation_rad) - floorplan.orbit_y_global * sin(rotation_rad);
        oy_rot = floorplan.orbit_x_global * sin(rotation_rad) + floorplan.orbit_y_global * cos(rotation_rad);
        floorplan.orbit_x_global = ox_rot; floorplan.orbit_y_global = oy_rot;
    end

    for i = 1:length(floorplan.element_positions)
        ep = floorplan.element_positions(i);
        xs = ep.pos_start(1) * cos(rotation_rad) - ep.pos_start(2) * sin(rotation_rad);
        ys = ep.pos_start(1) * sin(rotation_rad) + ep.pos_start(2) * cos(rotation_rad);
        xc = ep.pos_center(1) * cos(rotation_rad) - ep.pos_center(2) * sin(rotation_rad);
        yc = ep.pos_center(1) * sin(rotation_rad) + ep.pos_center(2) * cos(rotation_rad);
        xe = ep.pos_end(1) * cos(rotation_rad) - ep.pos_end(2) * sin(rotation_rad);
        ye = ep.pos_end(1) * sin(rotation_rad) + ep.pos_end(2) * cos(rotation_rad);
        floorplan.element_positions(i).pos_start = [xs, ys];
        floorplan.element_positions(i).pos_center = [xc, yc];
        floorplan.element_positions(i).pos_end = [xe, ye];
    end

    for i = 1:length(floorplan.elements)
        e = floorplan.elements(i);
        xc = e.x_center * cos(rotation_rad) - e.y_center * sin(rotation_rad);
        yc = e.x_center * sin(rotation_rad) + e.y_center * cos(rotation_rad);
        floorplan.elements(i).x_center = xc; floorplan.elements(i).y_center = yc;
        floorplan.elements(i).theta_center = e.theta_center + rotation_rad;
        xs = e.pos_start(1) * cos(rotation_rad) - e.pos_start(2) * sin(rotation_rad);
        ys = e.pos_start(1) * sin(rotation_rad) + e.pos_start(2) * cos(rotation_rad);
        floorplan.elements(i).pos_start = [xs, ys];
    end
end

floorplan.x_min = min(floorplan.x); floorplan.x_max = max(floorplan.x);
floorplan.y_min = min(floorplan.y); floorplan.y_max = max(floorplan.y);

if doPlot ~= 1
    result = struct();
    result.xmin = floorplan.x_min; result.xmax = floorplan.x_max;
    result.ymin = floorplan.y_min; result.ymax = floorplan.y_max;
    result.element_positions = floorplan.element_positions;
    if options.plot_orbit
        result.orbit_s = floorplan.orbit_s;
        result.orbit_state = floorplan.orbit_state;
    end
    return;
end

figure('Position', [100, 100, options.figure_size], 'Color', 'white');
plot(floorplan.x, floorplan.y, 'k-', 'LineWidth', options.line_width); hold on;

legend_handles = []; legend_texts = {}; element_types_shown = {};

for i = 1:length(floorplan.elements)
    element = floorplan.elements(i);
    x_center = element.x_center; y_center = element.y_center;
    element_type = element.type; element_theta = element.theta_center;

    switch element_type
        case 'quadrupole'
            fill_color = [.6, 0, 0]; element_display_name = 'Quadrupole';
            h = drawRectangle(x_center, y_center, max(element.length*options.rect_width_factor, 0.05), options.rect_length, element_theta, fill_color);
        case 'dipole'
            fill_color = [0.1, 0, 0.8]; element_display_name = 'Dipole';
            if element.h ~= 0
                angle = -element.length * element.h;
                radius = element.length / angle;
                rotated_pos_start = element.pos_start;
                rotated_start_theta = element.theta_center - angle/2;
                rotated_circle_center_x = rotated_pos_start(1) - radius * sin(rotated_start_theta);
                rotated_circle_center_y = rotated_pos_start(2) + radius * cos(rotated_start_theta);
                rotated_circle_center = [rotated_circle_center_x, rotated_circle_center_y];
                start_geo_angle = rotated_start_theta - pi/2;
                h = drawCurvedRectangle(rotated_circle_center, start_geo_angle, angle, radius, options.rect_length, fill_color);
            else
                rect_width = max(element.length*options.rect_width_factor, 0.05);
                h = drawRectangle(x_center, y_center, rect_width, options.rect_length, element_theta, fill_color);
            end
        case 'sextupole'
            fill_color = 'g'; element_display_name = 'Sextupole';
            h = drawRectangle(x_center, y_center, max(element.length*options.rect_width_factor, 0.05), options.rect_length, element_theta, fill_color);
        case 'cavity'
            fill_color = [0.5, 0, 0.5]; element_display_name = 'Cavity';
            num_cells = element.num_cells; cell_length = element.length / num_cells;
            h_group = hggroup;
            for j = 1:num_cells
                cell_offset = (j - 0.5) * cell_length;
                cx = element.pos_start(1) + cell_offset * cos(element_theta);
                cy = element.pos_start(2) + cell_offset * sin(element_theta);
                h_cell = drawEllipse(cx, cy, cell_length, options.rect_length * 1.2, element_theta, fill_color);
                set(h_cell, 'Parent', h_group);
            end
            h = h_group;
        case 'TGU'
            element_display_name = 'Undulator';
            num_segments = 9; segment_length = element.length / num_segments;
            colors = {[0.8, 0, 0], [0, 0, 0]};
            h_group = hggroup;
            for j = 1:num_segments
                segment_offset = (j - 0.5) * segment_length;
                cx = element.pos_start(1) + segment_offset * cos(element_theta);
                cy = element.pos_start(2) + segment_offset * sin(element_theta);
                h_segment = drawRectangle(cx, cy, segment_length*0.9, options.rect_length, element_theta, colors{mod(j-1,2)+1});
                set(h_segment, 'Parent', h_group);
            end
            h = h_group;
        otherwise
            % ====== 改动：未知非drift元件自动分配颜色 ======
            type_lower = lower(element_type);
            if auto_color_map.isKey(type_lower)
                fill_color = auto_color_map(type_lower);
            else
                auto_color_idx = auto_color_idx + 1;
                cidx = mod(auto_color_idx - 1, size(auto_color_palette, 1)) + 1;
                fill_color = auto_color_palette(cidx, :);
                auto_color_map(type_lower) = fill_color;
            end
            % 显示名称：首字母大写
            element_display_name = [upper(element_type(1)), lower(element_type(2:end))];
            h = drawRectangle(x_center, y_center, max(element.length*options.rect_width_factor, 0.05), options.rect_length, element_theta, fill_color);
    end

    if ~ismember(element_type, element_types_shown)
        element_types_shown{end+1} = element_type;
        legend_handles(end+1) = h; legend_texts{end+1} = element_display_name;
    end
end

if options.plot_orbit
    plot(floorplan.orbit_x_global, floorplan.orbit_y_global, 'Color', options.orbit_color, 'LineWidth', 1.5, 'LineStyle', '-');
    h_orbit = plot(NaN, NaN, 'Color', options.orbit_color, 'LineWidth', 1.5, 'LineStyle', '-');
    legend_handles(end+1) = h_orbit;
    legend_texts{end+1} = sprintf('Beam Orbit (x%d)', options.orbit_scale);
end

if ~isempty(legend_handles)
    lgd = legend(legend_handles, legend_texts, 'Location', 'northeast');
    set(lgd, 'FontName', options.font_name, 'FontSize', options.font_size);
end

title(options.title_text, 'FontSize', options.title_font_size, 'FontName', options.font_name);
xlabel('X (m)'); ylabel('Y (m)');
axis equal; box on;
set(gca, 'FontName', options.font_name, 'FontSize', options.font_size, 'LineWidth', options.axis_line_width);
xlim([floorplan.x_min - 2, floorplan.x_max + 2]);
ylim([floorplan.y_min - 2, floorplan.y_max + 2]);

result = floorplan;
end

function h = drawRectangle(x_center, y_center, rect_width, rect_length, theta, fill_color)
    perp_dir = theta + pi/2;
    p1_x = x_center + (rect_length/2) * cos(perp_dir) + (rect_width/2) * cos(theta);
    p1_y = y_center + (rect_length/2) * sin(perp_dir) + (rect_width/2) * sin(theta);
    p2_x = x_center + (rect_length/2) * cos(perp_dir) - (rect_width/2) * cos(theta);
    p2_y = y_center + (rect_length/2) * sin(perp_dir) - (rect_width/2) * sin(theta);
    p3_x = x_center - (rect_length/2) * cos(perp_dir) - (rect_width/2) * cos(theta);
    p3_y = y_center - (rect_length/2) * sin(perp_dir) - (rect_width/2) * sin(theta);
    p4_x = x_center - (rect_length/2) * cos(perp_dir) + (rect_width/2) * cos(theta);
    p4_y = y_center - (rect_length/2) * sin(perp_dir) + (rect_width/2) * sin(theta);
    h = patch([p1_x, p2_x, p3_x, p4_x], [p1_y, p2_y, p3_y, p4_y], fill_color, 'EdgeColor', 'none');
end

function h = drawCurvedRectangle(circle_center, start_angle, angle, R, thickness, fill_color)
    num_samples = 30;
    end_angle = start_angle + angle;
    angles = linspace(start_angle, end_angle, num_samples);
    r_outer = R + thickness/2; r_inner = R - thickness/2;
    x_out = circle_center(1) + r_outer * cos(angles);
    y_out = circle_center(2) + r_outer * sin(angles);
    x_in = circle_center(1) + r_inner * cos(flip(angles));
    y_in = circle_center(2) + r_inner * sin(flip(angles));
    x_global = [x_out, x_in]; y_global = [y_out, y_in];
    h = patch(x_global, y_global, fill_color, 'EdgeColor', 'none');
end

function h = drawEllipse(x_center, y_center, width, height, theta, fill_color)
    t = linspace(0, 2*pi, 40);
    x_local = (width/2) * cos(t); y_local = (height/2) * sin(t);
    x_global = x_center + x_local * cos(theta) - y_local * sin(theta);
    y_global = y_center + x_local * sin(theta) + y_local * cos(theta);
    h = patch(x_global, y_global, fill_color, 'EdgeColor', 'none');
end