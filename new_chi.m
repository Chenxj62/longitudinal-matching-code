function [y, cons] = new_chi(x)
rho = 2;
cons = [];  % 默认初始化约束值

Lq = .2;

% % 设置高精度模式
% format long;
% 
% 使用额外精度的系数定义
x_high = x;


% 定义四极铁元件
ELE.Q1 = createElement('quadrupole', Lq, x_high(1), 'Q1');
ELE.Q2 = createElement('quadrupole', Lq, x_high(2), 'Q2');
ELE.Q3 = createElement('quadrupole', Lq, x_high(3), 'Q3');
ELE.Q4 = createElement('quadrupole', Lq, x_high(4), 'Q4');
ELE.Q5 = createElement('quadrupole', Lq/2, x_high(5), 'Q5');


% 定义漂移段
ELE.D1 = createElement('drift', x(6), 0, 'D1');
ELE.D2 = createElement('drift', x(7), 0, 'D2');
ELE.D3 = createElement('drift', x(8), 0, 'D3');
ELE.D4 = createElement('drift', x(9), 0, 'D4');
ELE.D5 = createElement('drift', x(10), 0, 'D5');
ELE.D6 = createElement('drift', x(11), 0, 'D6');

k1=x(12);
rho2=k1*rho;

ELE.S1 = createElement('sextupole', 0.1, x(13), 'S1');
ELE.S2 = createElement('sextupole', 0.1, x(14), 'S2');
ELE.S3 = createElement('sextupole', 0.1, x(15), 'S3');







% 定义偏转铁
ELE.B1 = createElement('dipole', x(16)*rho, 1/rho, 'B1');
ELE.B2 = createElement('dipole', x(16)*rho2, -1/rho2, 'B2');
% 定义偏转铁

% 定义每段磁路 - 每段4个四级铁，首尾都是漂移段
latm1 = [ELE.D1, ELE.Q1, ELE.D2,ELE.S1,ELE.D2 ELE.Q2, ELE.D3,ELE.S2,ELE.D3 ELE.Q3, ELE.D4];

latm2 = [ELE.D5, ELE.Q4, ELE.D6,ELE.S3,ELE.D6, ELE.Q5];

midbeamline=[ELE.B1,latm1,ELE.B2,latm2];

neglatm1=reverseBeamline(midbeamline);




T266=nonlinear_calculator(midbeamline,2,6,6,0.0001);




beamline = [midbeamline,neglatm1];

[R11min,R11max,~,~]=findRijmax(beamline,1,1);
[R12min,R12max,~,~]=findRijmax(beamline,1,2);
[R21min,R21max,~,~]=findRijmax(beamline,2,1);
[R22min,R22max,~,~]=findRijmax(beamline,2,2);
[R33min,R33max,~,~]=findRijmax(beamline,3,3);
[R34min,R34max,~,~]=findRijmax(beamline,3,4);
[R43min,R43max,~,~]=findRijmax(beamline,4,3);
[R44min,R44max,~,~]=findRijmax(beamline,4,4);

Rmax = max([abs(R11max), abs(R11min), ...
            abs(R12max), abs(R12min), ...
            abs(R21max), abs(R21min), ...
            abs(R22max), abs(R22min), ...
            abs(R33max), abs(R33min), ...
            abs(R34max), abs(R34min), ...
            abs(R43max), abs(R43min), ...
            abs(R44max), abs(R44min)]);

T566=nonlinear_calculator(midbeamline,5,6,6,0.0001);

Tnon=0;%abs(T266);%max(abs(T266),abs(T566));
r1=calculateTotalMatrix(midbeamline);  %%%得到一半的传输矩阵

r=calculateTotalMatrix(beamline);
% 其余代码不变

[intB,intD,~]=CalcCSRkick(beamline,30,0,1e-4);




if abs(r1(2,6))<1e-4 
    cdis=1;
else
    cdis=(abs(r1(2,6))/0.0001);
end
f1=abs((intB(1)))+abs((intB(2)));
f2=abs(intD(1))+abs(intD(2));



 if Tnon<1E-2

c566=1;
 else
     c566=abs(Tnon/1e-2);
 end

if abs(r(1,1))<3&&abs(r(3,3))<3
c16=1;
else
    c16=((abs(r(1,1))+abs(r(3,3)))/0.1);
end

if Rmax<15
    c1 = 1;
else
    c1 = (abs(Rmax)/.1);
end

if r(5,6)<-0.4||r(5,6)>-0.15
c2=abs(r(5,6)/1e-3)+2;
else
    c2=1;
end

y = [c16*cdis*c1*c2*c566*f1+c16+cdis+c1+c2+c566-5, c16*cdis*c1*c2*c566*f2+c16+cdis+c1+c2+c566-5];


end