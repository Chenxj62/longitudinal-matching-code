
function reversed_beamline = rev_beamline(beamline)
    % 初始化空数组
    reversed_beamline = [];
    
    % 反向排列元件
    for i = length(beamline):-1:1

        if strcmpi(beamline(i).type, 'cavity')

            phi = beamline(i).phase+180;

            beamline(i).phase=phi;

        end
        % 复制原元件
        element = beamline(i);
        reversed_element = element;
        
        % 将元件添加到反向束线
        reversed_beamline = [reversed_beamline, reversed_element];
    end
    
    % 注意：这里保持所有元件的强度符号不变
    % 如果需要改变某些元件的强度符号，可以在这里添加代码
end