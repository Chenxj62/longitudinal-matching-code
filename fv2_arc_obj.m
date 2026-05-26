function [y, cons] = fv2_arc_obj(x)



cons = [];                 % 约束初始化


adjust_vars = x;

% 基础参数
B1_angle = 10*pi/180;  % rad
rho1_base = 3.45721;  % B1基础弯转半径
B2_angle = 10*pi/180;  % rad
rho2_base = 4.26259;  % B2基础弯转半径
B3_angle = 10*pi/180;  % rad
rho3_base = 4.97707;  % B3基础弯转半径
B4_angle = 10*pi/180;  % rad
rho4_base = 4.58811;  % B4基础弯转半径

% 四极铁长度
Lq1 = 0.300000;  % Q1长度
Lq2 = 0.300000;  % Q2长度
Lq3 = 0.300000;  % Q3长度
Lq4 = 0.300000;  % Q4长度
Lq5 = 0.300000;  % Q5长度
Lq6 = 0.300000;  % Q6长度
Lq7 = 0.300000;  % Q7长度

cons = [];  % 约束初始化
eang = 5*pi/180;

% 计算弯铁投影长度（chord length），精确到毫米
% Chord length = 2*rho*sin(theta/2)
chord_B1 = round(2 * (rho2_base + adjust_vars(20)) * sin(B1_angle/2), 2);
chord_B2 = round(2 * (rho2_base + adjust_vars(20)) * sin(B2_angle/2), 2);
chord_B3 = round(2 * (rho2_base + adjust_vars(20)) * sin(B3_angle/2), 2);
chord_B4 = round(2 * (rho2_base + adjust_vars(20)) * sin(B4_angle/2), 2);

% 从投影长度反推实际弯转半径
rho1_actual = chord_B1 / (2 * sin(B1_angle/2));
rho2_actual = chord_B2 / (2 * sin(B2_angle/2));
rho3_actual = chord_B3 / (2 * sin(B3_angle/2));
rho4_actual = chord_B4 / (2 * sin(B4_angle/2));

% 计算弧长（dipole长度）
L_B1 = rho1_actual * B1_angle;
L_B2 = rho2_actual * B2_angle;
L_B3 = rho3_actual * B3_angle;
L_B4 = rho4_actual * B4_angle;

% 创建元件
ELE.B1 = createElement('dipole', L_B1, B1_angle/L_B1, 'B1');
E1 = createElement('edge', 0, B1_angle/L_B1, eang, ' ');

ELE.D1 = createElement('drift', round(0.514706 + adjust_vars(8), 3), 0, 'D1');
ELE.Q1 = createElement('quadrupole', Lq1, 3.89957 + adjust_vars(1), 'QM1');

ELE.D2 = createElement('drift', round(0.815087 + adjust_vars(9), 3), 0, 'D2');
ELE.Q2 = createElement('quadrupole', Lq2, -3.90221 + adjust_vars(2), 'QM2');

ELE.D3 = createElement('drift', round(1.509774 + adjust_vars(10), 3), 0, 'D3');
ELE.Q3 = createElement('quadrupole', Lq3, 3.35649 + adjust_vars(3), 'QM3');

ELE.D4 = createElement('drift', round(0.383658 + adjust_vars(11), 3), 0, 'D4');
ELE.B2 = createElement('dipole', L_B2, B2_angle/L_B2, 'B2');
E2 = createElement('edge', 0, B2_angle/L_B2, eang, ' ');

ELE.D5 = createElement('drift', round(0.564785 + adjust_vars(12), 3), 0, 'D5');
ELE.Q4 = createElement('quadrupole', Lq4, -3.36332 + adjust_vars(4), 'Q4');

ELE.D6 = createElement('drift', round(1.080429 + adjust_vars(13), 3), 0, 'D6');
ELE.D7 = createElement('drift', round(1.181531 + adjust_vars(14), 3), 0, 'D7');
ELE.Q5 = createElement('quadrupole', Lq5, 2.27378 + adjust_vars(5), 'Q5');

ELE.D8 = createElement('drift', round(0.809512 + adjust_vars(15), 3), 0, 'D8');
ELE.B3 = createElement('dipole', L_B3, B3_angle/L_B3, 'B3');
E3 = createElement('edge', 0, B3_angle/L_B3, eang, ' ');

ELE.D9 = createElement('drift', round(0.967255 + adjust_vars(16), 3), 0, 'D9');
ELE.Q6 = createElement('quadrupole', Lq6, -2.60891 + adjust_vars(6), 'Q6');

ELE.D10 = createElement('drift', round(0.240294 + adjust_vars(17), 3), 0, 'D10');
ELE.D11 = createElement('drift', round(0.964294 + adjust_vars(18), 3), 0, 'D11');
ELE.Q7 = createElement('quadrupole', Lq7, 4.12778 + adjust_vars(7), 'Q7');

ELE.D12 = createElement('drift', round(0.741644 + adjust_vars(19), 3), 0, 'D12');
ELE.B4 = createElement('dipole', L_B4, B4_angle/L_B4, 'B4');
E4 = createElement('edge', 0, B4_angle/L_B4, eang, ' ');
% 定义完整束线
beamline = [E1, ELE.B1, E1, ELE.D1, ELE.Q1, ELE.D2, ELE.Q2, ...
                     ELE.D3, ELE.Q3, ELE.D4, E2, ELE.B2, E2, ELE.D5, ...
                     ELE.Q4, ELE.D6, ELE.D7, ELE.Q5, ELE.D8, ...
                     E3, ELE.B3, E3, ELE.D9, ELE.Q6, ELE.D10, ELE.D11, ...
                     ELE.Q7, ELE.D12, E4, ELE.B4, E4];







beamline1 = reverseBeamline(beamline);



[R11min,R11max,~,~]=findRijmax(beamline,1,1);
[R12min,R12max,~,~]=findRijmax(beamline,1,2);
[R21min,R21max,~,~]=findRijmax(beamline,2,1);
[R22min,R22max,~,~]=findRijmax(beamline,2,2);
[R33min,R33max,~,~]=findRijmax(beamline,3,3);
[R34min,R34max,~,~]=findRijmax(beamline,3,4);
[R43min,R43max,~,~]=findRijmax(beamline,4,3);
[R44min,R44max,~,~]=findRijmax(beamline,4,4);

Rmax = max([abs(R33max), abs(R33min), ...
            abs(R34max), abs(R34min), ...
            abs(R43max), abs(R43min), ...
            abs(R44max), abs(R44min)]);

Rmax1 = max([abs(R11max), abs(R11min), ...
            abs(R12max), abs(R12min), ...
            abs(R21max), abs(R21min), ...
            abs(R22max), abs(R22min)]);



r=calculateTotalMatrix(beamline);
% 其余代码不变

[intB,intD,~]=CalcCSRkick(beamline1,250,0,1.0e-5);




if abs(r(2,6))<1e-3 && abs(r(1,6))<1e-3
    cdis=1;
else
    cdis=(abs(r(2,6))/1e-3)+(abs(r(1,6))/1e-3);
end




 if r(5,6)<-.25

  c566=exp(abs(r(5,6)/2e-1));

 elseif r(5,6)>-0.07
     c566=exp(abs(r(5,6)+1));
 else
c566=1;
 end


c16=1;



f1=abs(intB(1));%+abs((intB(2)));
f2=abs(intB(2));%+abs((intD(2)));

if Rmax<5
    c1 = 1;
else
    c1 = (abs(Rmax)/.1);
end


if Rmax1<5
    c2 = 1;
else
    c2 = (abs(Rmax1)/.1);
end

y = [c16*cdis*c1*c2*c566*f1+c16+cdis+c1+c2+c566-5, c16*cdis*c1*c2*c566*f2+c16+cdis+c1+c2+c566-5];


end