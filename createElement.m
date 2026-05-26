function element = createElement(type, length, varargin)
    element.type = type;
    element.length = length;
    
    % 默认值
    element.k1 = 0;     % 四极场强度 (m^-2)
    element.k2 = 0;     % 六极场强度 (m^-3)
    element.h = 0;      % 曲率 = 1/rho (m^-1)
    element.angle = 0;  % 偏转角度 (rad)
    element.name = '';  % 元件名称
    
    % Cavity特有参数
    element.gradient = 0;     % 腔梯度 (MV/m)
    element.num_cells = 0;    % cell数量
    element.frequency = 0;    % 频率 (MHz)
    element.phase = 0;        % 相位 (度)

    %TGU 特有参数
    element.gamma=0; %能量
    element.alp=0;  %磁场梯度
    element.alpc=0; %矫正线圈梯度
    element.k0=0;
    element.ku=0;

    % Matrix元件参数 (自定义传输矩阵)
    % 辛条件: m11*m22 - m12*m21 = 1 => m22 = (1 + m12*m21) / m11
    %         m33*m44 - m34*m43 = 1 => m44 = (1 + m34*m43) / m33
    % 纵向: 单位矩阵
    element.m11 = 1;
    element.m22 = 1;
    element.m21 = 0;
    element.m33 = 1;
    element.m44 = 1;
    element.m43 = 0;
    
    %%1st order Taylor Map
    

    % 根据元件类型设置参数
    if ~isempty(varargin)
        switch type
            case 'edge'
               
                element.h=varargin{1};
                element.angle=varargin{2};
            case 'dipole'
                element.h = varargin{1};  % 曲率
                if abs(element.h) > 1e-10
                    element.angle = length * element.h;
                end
                % 检查是否有名称参数
                if numel(varargin) >= 2
                    element.name = varargin{2};
                end
                
            case {'quadrupole','quad'}
                element.k1 = varargin{1};  % 四极场强度
                % 检查是否有名称参数
                if numel(varargin) >= 2
                    element.name = varargin{2};
                end
                
            case {'sextupole','sext'}
                element.k2 = varargin{1};  % 六极场强度
                % 检查是否有名称参数
                if numel(varargin) >= 2
                    element.name = varargin{2};
                end
                
            case 'cavity'
                % Cavity参数处理，更加健壮的方式
                if numel(varargin) >= 1
                    if isnumeric(varargin{1}) && numel(varargin{1}) >= 4
                        % 数组形式参数: [梯度, cell数, 频率, 相位]
                        params = varargin{1};
                        element.gradient = params(1);     % 梯度 (MV/m)
                        element.num_cells = params(2);    % cell数量
                        element.frequency = params(3);    % 频率 (MHz)
                        element.phase = params(4);        % 相位 (度)
                          element.length=element.num_cells*0.5*3e8/(element.frequency*1e6);
                        % 检查是否有名称参数
                        if numel(varargin) >= 2
                            element.name = varargin{2};
                        end
                    elseif isnumeric(varargin{1})
                        % 单独参数形式，检查是否有足够的参数
                        if numel(varargin) >= 4
                            element.gradient = varargin{1};    % 梯度 (MV/m)
                            element.num_cells = varargin{2};   % cell数量
                            element.frequency = varargin{3};   % 频率 (MHz)
                            element.phase = varargin{4};       % 相位 (度)
                            
                            % 检查是否有名称参数
                            if numel(varargin) >= 5
                                element.name = varargin{5};
                            end
                        else
                            warning('Cavity needs at least 4 parameters: gradient, num_cells, frequency, phase');
                        end
                    end
                end
                
            case 'marker'
                if ~isempty(varargin)
                    element.name = varargin{1};
                end

            case '1stTTay'
                % 自定义传输矩阵元件
                % 用法1: createElement('matrix', L, [m11,m12,m21,m33,m34,m43], 'name')
                
                %
                % 辛条件自动求解 m22, m44:
                %   m22 = (1 + m12*m21) / m11
                %   m44 = (1 + m34*m43) / m33
                % 纵向 (5,6) 平面: 单位矩阵
                
                if numel(varargin) >= 1
                    if isnumeric(varargin{1}) && numel(varargin{1}) >= 6
                        % 数组形式: [m11, m12, m21, m33, m34, m43]
                        params = varargin{1};
                        element.m11 = params(1);
                        element.m22 = params(2);
                        element.m21 = params(3);
                        element.m33 = params(4);
                        element.m44 = params(5);
                        element.m43 = params(6);

                        eps_m = 1e-12;

                        if abs(element.m21) < eps_m
                            if element.m21 < 0
                                element.m21 = -eps_m;
                            else
                                element.m21 = eps_m;
                            end
                        end
                        
                        if abs(element.m43) < eps_m
                            if element.m43 < 0
                                element.m43 = -eps_m;
                            else
                                element.m43 = eps_m;
                            end
                        end
                        
                        if numel(varargin) >= 2
                            element.name = varargin{2};
                        end
                    elseif isnumeric(varargin{1}) && numel(varargin) >= 6
                        % 单独参数形式
                        element.m11 = varargin{1};
                        element.m22 = varargin{2};
                        element.m21 = varargin{3};
                        element.m33 = varargin{4};
                        element.m44 = varargin{5};
                        element.m43 = varargin{6};
                        
                        if numel(varargin) >= 7
                            element.name = varargin{7};
                        end
                    else
                        warning('Matrix element needs 6 parameters: m11, m12, m21, m33, m34, m43');
                    end
                end
                
                % 检查辛条件可解性，不可解则回退为单位矩阵
                if abs(element.m11) < 1e-15 || abs(element.m33) < 1e-15
                    if abs(element.m11) < 1e-15
                        warning('matrix element: m11 ≈ 0, symplecticity unsolvable. Falling back to identity.');
                         element.m11 = 1; element.m12 = 0; element.m21 = 0;
                    end
                    if abs(element.m33) < 1e-15
                        warning('matrix element: m33 ≈ 0, symplecticity unsolvable. Falling back to identity.');
                        element.m33 = 1; element.m34 = 0; element.m43 = 0;
                    end
                   
                    
                end

            case 'TGU'
                element.gamma=varargin{1};
                 element.alp=varargin{2};  %磁场梯度
                 element.alpc=varargin{3}; %矫正线圈梯度
                 element.k0=varargin{4};
                element.ku=varargin{5};
                element.name=varargin{6};
            case 'drift'
                % 对漂移段的可选名称参数
                if ~isempty(varargin)
                    element.name = varargin{2};
                end
        end
    end
end