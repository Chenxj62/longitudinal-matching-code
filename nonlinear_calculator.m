function varargout = nonlinear_calculator(beamline, varargin)
    % 计算非线性传输矩阵元素 Tijk，支持两种调用方式:
    %
    % 【单点模式】返回束线末端的一个 Tijk 值:
    %   Tijk = nonlinear_calculator(beamline, i, j, k)
    %   Tijk = nonlinear_calculator(beamline, i, j, k, delta)
    %   Tijk = nonlinear_calculator(beamline, i, j, k, delta, E0)
    %
    % 【演化模式】返回指定 Tijk 项沿束线的演化 (O(n)单次遍历):
    %   [s, T] = nonlinear_calculator(beamline, ilist, jlist, klist)
    %   [s, T] = nonlinear_calculator(beamline, ilist, jlist, klist, delta)
    %   [s, T] = nonlinear_calculator(beamline, ilist, jlist, klist, delta, E0)
    %
    %   s - 纵向位置数组 (n+1)×1
    %   T - 演化矩阵 (n+1)×m, 每列对应 ilist/jlist/klist 中的一个 Tijk 项
    %
    % 公共参数:
    %   delta - 计算偏差值 (默认 0.001)
    %   E0    - 初始能量 eV  (默认 4e6)
    %
    % 示例:
    %   T566 = nonlinear_calculator(beamline, 5, 6, 6, 0.01);
    %
    %   [s, T] = nonlinear_calculator(beamline, [5,1], [6,6], [6,6], 0.001, E0);
    %   T566_evol = T(:,1);
    %   T166_evol = T(:,2);

    delta = 0.001;
    E0    = 4e6;
    if nargin >= 6, E0    = varargin{5}; end
    if nargin >= 5, delta = varargin{4}; end

    ilist = varargin{1}(:)';
    jlist = varargin{2}(:)';
    klist = varargin{3}(:)';

    if nargout <= 1
        % ── 单点模式 ────────────────────────────────────────────────
        varargout{1} = single_Tijk(beamline, ilist(1), jlist(1), klist(1), delta, E0);
    else
        % ── 演化模式 ────────────────────────────────────────────────
        [varargout{1}, varargout{2}] = evolution(beamline, ilist, jlist, klist, delta, E0);
    end
end

% =========================================================================
% 单点模式: 计算束线末端的 Tijk
% =========================================================================
function Tijk = single_Tijk(beamline, iindex, jindex, kindex, delta, E0)
    if jindex == kindex
        % 对角项: 3粒子公式 [f(+δ) + f(-δ) - 2f(0)] / (2δ²)
        p0 = zeros(6,1);  p1 = zeros(6,1);  p2 = zeros(6,1);
        p1(jindex) =  delta;
        p2(jindex) = -delta;
        p0 = track_beamline(p0, beamline, E0);
        p1 = track_beamline(p1, beamline, E0);
        p2 = track_beamline(p2, beamline, E0);
        Tijk = (p1(iindex) + p2(iindex) - 2*p0(iindex)) / (2*delta^2);
    else
        % 非对角项: 4粒子公式 [f(++) - f(+-) - f(-+) + f(--)] / (4δ²)
        p1 = zeros(6,1);  p2 = zeros(6,1);
        p3 = zeros(6,1);  p4 = zeros(6,1);
        p1(jindex) =  delta;  p1(kindex) =  delta;
        p2(jindex) =  delta;  p2(kindex) = -delta;
        p3(jindex) = -delta;  p3(kindex) =  delta;
        p4(jindex) = -delta;  p4(kindex) = -delta;
        p1 = track_beamline(p1, beamline, E0);
        p2 = track_beamline(p2, beamline, E0);
        p3 = track_beamline(p3, beamline, E0);
        p4 = track_beamline(p4, beamline, E0);
        Tijk = (p1(iindex) - p2(iindex) - p3(iindex) + p4(iindex)) / (4*delta^2);
    end
end

% =========================================================================
% 演化模式: O(n)单次遍历，按(j,k)分组追踪粒子
% =========================================================================
function [s_positions, T_evol] = evolution(beamline, ilist, jlist, klist, delta, E0)
    m = length(ilist);

    % 找出唯一(j,k)对，将请求的项映射到组
    jk_pairs = [jlist(:), klist(:)];
    [unique_jk, ~, term2group] = unique(jk_pairs, 'rows');
    num_groups = size(unique_jk, 1);

    % 为每个唯一(j,k)对初始化粒子矩阵 (6×Np)
    group_p = cell(num_groups, 1);
    for g = 1:num_groups
        j = unique_jk(g, 1);
        k = unique_jk(g, 2);
        if j == k
            % 对角项: 3粒子 [p0 | p+ | p-]
            p0 = zeros(6,1);
            pp = zeros(6,1);  pp(j) =  delta;
            pm = zeros(6,1);  pm(j) = -delta;
            group_p{g} = [p0, pp, pm];
        else
            % 非对角项: 4粒子 [p++ | p+- | p-+ | p--]
            p1 = zeros(6,1);  p1(j)= delta;  p1(k)= delta;
            p2 = zeros(6,1);  p2(j)= delta;  p2(k)=-delta;
            p3 = zeros(6,1);  p3(j)=-delta;  p3(k)= delta;
            p4 = zeros(6,1);  p4(j)=-delta;  p4(k)=-delta;
            group_p{g} = [p1, p2, p3, p4];
        end
    end

    n = length(beamline);
    s_positions = zeros(n+1, 1);
    T_evol      = zeros(n+1, m);
    current_s   = 0;
    current_E   = E0;

    for i = 1:n
        element = beamline(i);

        % 追踪每组所有粒子穿过当前元件
        for g = 1:num_groups
            np = size(group_p{g}, 2);
            for p = 1:np
                group_p{g}(:,p) = track_element(group_p{g}(:,p), element, current_E);
            end
        end

        if strcmpi(element.type, 'cavity')
            current_E = update_cavity_energy(element, current_E);
        end

        % 计算所有请求的 Tijk
        for t = 1:m
            g  = term2group(t);
            ii = ilist(t);
            j  = unique_jk(g, 1);
            k  = unique_jk(g, 2);
            P  = group_p{g};
            if j == k
                T_evol(i+1, t) = (P(ii,2) + P(ii,3) - 2*P(ii,1)) / (2*delta^2);
            else
                T_evol(i+1, t) = (P(ii,1) - P(ii,2) - P(ii,3) + P(ii,4)) / (4*delta^2);
            end
        end

        if isfield(element, 'length')
            current_s = current_s + element.length;
        end
        s_positions(i+1) = current_s;
    end
end

% =========================================================================
% 追踪粒子穿过整段束线
% =========================================================================
function p = track_beamline(p, beamline, E0)
    current_E = E0;
    for i = 1:length(beamline)
        p = track_element(p, beamline(i), current_E);
        if strcmpi(beamline(i).type, 'cavity')
            current_E = update_cavity_energy(beamline(i), current_E);
        end
    end
end

% =========================================================================
% 追踪粒子穿过单个元件
% =========================================================================
function p = track_element(p, element, current_E)
    me_c2 = 0.511e6;

    if strcmpi(element.type, 'drift')
        p = trackDrift(p, element.length);

    elseif strcmpi(element.type, 'marker')
        % 无操作

    elseif strcmpi(element.type, 'edge')
        p = trackEdge(p, element.h, element.angle);

    elseif strcmpi(element.type, 'quadrupole')
        if abs(element.k1) < 1e-4
            p = trackDrift(p, element.length);
        else
            p = trackQuadrupole(p, element.length, element.k1);
        end

    elseif strcmpi(element.type, 'sextupole')
        p = trackSextupole(p, element.length, element.k2);

    elseif strcmpi(element.type, 'dipole')
        if isfield(element, 'h') && abs(element.h) > 1e-10
            g = element.h;
        else
            g = 1/5.0;
            warning('找不到有效的曲率字段h，使用默认半径5.0米');
        end
        if isfield(element, 'angle') && element.angle ~= 0
            angle = element.angle;
        else
            angle = g * element.length;
        end
        k1 = 0;
        if isfield(element, 'k1'), k1 = element.k1; end
        p = trackDipole(p, g, (1/g)*angle, k1);

    elseif strcmpi(element.type, '1stTTay')
        R = getElementMatrix(element);
        p = R * p;

    elseif strcmpi(element.type, 'cavity')
        gamma_in = current_E / me_c2;
        R = getElementMatrix(element, gamma_in);
        p = R * p;

    else
        % 未识别元件 (kicker / monitor / instrument / ecollimator / ...):
        % 按漂移段处理, 避免丢掉长度引起的横向位移误差。
        if isfield(element, 'length') && element.length > 0
            p = trackDrift(p, element.length);
        end
    end
end

function new_E = update_cavity_energy(element, current_E)
    if isfield(element, 'num_cells') && isfield(element, 'frequency')
        L_cav = element.num_cells * 0.5 * 3e8 / (element.frequency * 1e6);
    else
        L_cav = element.length;
    end
    A   = 0;  if isfield(element, 'gradient'), A   = element.gradient * 1e6;  end
    phi = 0;  if isfield(element, 'phase'),    phi = element.phase * pi/180;  end
    new_E = current_E + A * L_cav * cos(phi);
end

% =========================================================================
% 元件跟踪子函数
% =========================================================================
function p_final = trackDrift(p_initial, L)
    x1=p_initial(1); px1=p_initial(2); y1=p_initial(3);
    py1=p_initial(4); z1=p_initial(5); d=p_initial(6);
    pl = sqrt(1 - (px1^2 + py1^2)/(1+d)^2);
    p_final = [x1 + L*px1/((1+d)*pl);
               px1;
               y1 + L*py1/((1+d)*pl);
               py1;
               z1 + (1 - 1/pl)*L;
               d];
end

function p_final = trackEdge(p_initial, h, angle)
    p_final    = p_initial;
    p_final(2) = p_initial(2) + p_initial(1)*h*tan(angle);
    p_final(4) = p_initial(4) - p_initial(3)*h*tan(angle);
end

function p_final = trackQuadrupole(p_initial, L, k1)
    x1=p_initial(1); px1=p_initial(2); y1=p_initial(3);
    py1=p_initial(4); z1=p_initial(5); d=p_initial(6);
    wx = sqrt(abs(k1)/(1+d));  wy = wx;
    if k1 > 0
        cx=cos(wx*L);  sx=sin(wx*L)/wx;  tx=-1;
        cy=cosh(wy*L); sy=sinh(wy*L)/wy; ty= 1;
    else
        cx=cosh(wx*L); sx=sinh(wx*L)/wx; tx= 1;
        cy=cos(wy*L);  sy=sin(wy*L)/wy;  ty=-1;
    end
    x2  = cx*x1 + sx*px1/(1+d);
    px2 = tx*wx^2*(1+d)*sx*x1 + cx*px1;
    y2  = cy*y1 + sy*py1/(1+d);
    py2 = ty*wy^2*(1+d)*sy*y1 + cy*py1;
    z2  = z1 + (tx*wx^2*(L-cx*sx)/4)*x1^2 ...
             + (-tx*wx^2*sx^2/(2*(1+d)))*x1*px1 ...
             + (-(L+cx*sx)/(4*(1+d)^2))*px1^2 ...
             + (ty*wy^2*(L-cy*sy)/4)*y1^2 ...
             + (-ty*wy^2*sy^2/(2*(1+d)))*y1*py1 ...
             + (-(L+cy*sy)/(4*(1+d)^2))*py1^2;
    p_final = [x2; px2; y2; py2; z2; d];
end

function p_final = trackSextupole(p_initial, L, k2)
    x=p_initial(1); px=p_initial(2); y=p_initial(3);
    py=p_initial(4); z=p_initial(5); d=p_initial(6);
    p_final = p_initial;
    Lh = L/2;
    [x,y,z] = driftStep(x,px,y,py,z,d,Lh);
    p_final(1)=x; p_final(3)=y; p_final(5)=z;
    p_final(2) = p_final(2) - 0.5*k2*(p_final(1)^2 - p_final(3)^2)*L;
    p_final(4) = p_final(4) + k2*p_final(1)*p_final(3)*L;
    [x,y,z] = driftStep(p_final(1),p_final(2),p_final(3),p_final(4),p_final(5),d,Lh);
    p_final(1)=x; p_final(3)=y; p_final(5)=z;
end

function [xo,yo,zo] = driftStep(x,px,y,py,z,d,L)
    pl = sqrt(1-(px^2+py^2)/(1+d)^2);
    xo = x + L*px/((1+d)*pl);
    yo = y + L*py/((1+d)*pl);
    zo = z + (1-1/pl)*L;
end

function p_final = trackDipole(p_initial, g, L, ~)
    % k1=0 弯铁
    x=p_initial(1); px=p_initial(2); y=p_initial(3);
    py=p_initial(4); z=p_initial(5); d=p_initial(6);
    kx = g^2;
    xc = (g*(1+d) - g) / kx;   % = d/g
    wx = sqrt(kx/(1+d));
    cx=cos(wx*L); sx=sin(wx*L)/wx; tx=-1;
    wy=1e-5; cy=1; sy=L; ty=1;
    m5   = -g*xc*L;
    m51  = -g*sx;
    m52  =  tx*g*(1-cx)/((1+d)*wx^2);
    m511 =  tx*wx^2*(L-cx*sx)/4;
    m512 = -tx*wx^2*sx^2/(2*(1+d));
    m522 = -(L+cx*sx)/(4*(1+d)^2);
    m533 =  ty*wy^2*(L-cy*sy)/4;
    m534 = -ty*wy^2*sy^2/(2*(1+d));
    m544 = -(L+cy*sy)/(4*(1+d)^2);
    dx = x - xc;
    p_final = [cx*dx + sx*px/(1+d) + xc;
               tx*wx^2*(1+d)*sx*dx + cx*px;
               cy*y + sy*py/(1+d);
               ty*wy^2*(1+d)*sy*y + cy*py;
               z+m5+m51*dx+m52*px+m511*dx^2+m512*dx*px+m522*px^2+m533*y^2+m534*y*py+m544*py^2;
               d];
end
