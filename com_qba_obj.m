function [y, cons] = com_qba_obj(x)
rho = 1;
cons = [];  % 默认初始化约束值




Lq = .2;
ang = pi/18;




ELE.Q1 = createElement('quadrupole', Lq, x(1), 'Q1');
ELE.Q2 = createElement('quadrupole', Lq, x(2), 'Q2');
ELE.Q3 = createElement('quadrupole', Lq, x(3), 'Q3');
ELE.Q4 = createElement('quadrupole', Lq, x(4), 'Q4');
ELE.Q5 = createElement('quadrupole', Lq/2, x(5), 'Q5');


ELE.D1 = createElement('drift', x(6),  0, 'D1');
ELE.D2 = createElement('drift', x(7),  0, 'D2');
ELE.D3 = createElement('drift', x(8),  0, 'D3');
ELE.D4 = createElement('drift', x(9),  0, 'D4');
ELE.D5 = createElement('drift', x(10), 0, 'D5');
ELE.D6 = createElement('drift', x(11), 0, 'D6');



% 定义标记点和偏转铁
E1     =createElement('edge', 0, 1/rho,ang/2, ' ');
E2     =createElement('edge', 0, -1/rho,-ang/2, ' ');
ELE.B1 = createElement('dipole', ang*rho, 1/rho, 'B1');
ELE.B2 = createElement('dipole', ang*rho, -1/rho, 'B2');
% 定义子束线
latm1 = [ELE.D1, ELE.Q1, ELE.D2, ELE.Q2, ELE.D3, ELE.Q3, ELE.D4];

latm2 = [ELE.D5, ELE.Q4, ELE.D6, ELE.Q5];


midbeamline = [E1, ELE.B1, E1,  latm1, E2,ELE.B2, E2, latm2 ];

    % BEAMLINE_END
   rebeamline=reverseBeamline(midbeamline);

   beamline=[midbeamline,rebeamline];

   rmid=calculateTotalMatrix(midbeamline);

   r=calculateTotalMatrix(beamline);

% [R11min,R11max,~,~]=findRijmax(beamline,1,1);
% [R12min,R12max,~,~]=findRijmax(beamline,1,2);
% [R21min,R21max,~,~]=findRijmax(beamline,2,1);
% [R22min,R22max,~,~]=findRijmax(beamline,2,2);
% [R33min,R33max,~,~]=findRijmax(beamline,3,3);
% [R34min,R34max,~,~]=findRijmax(beamline,3,4);
% [R43min,R43max,~,~]=findRijmax(beamline,4,3);
% [R44min,R44max,~,~]=findRijmax(beamline,4,4);
% 
% Rmax = max([abs(R11max), abs(R11min), ...
%             abs(R12max), abs(R12min), ...
%             abs(R21max), abs(R21min), ...
%             abs(R22max), abs(R22min), ...
%             abs(R33max), abs(R33min), ...
%             abs(R34max), abs(R34min), ...
%             abs(R43max), abs(R43min), ...
%             abs(R44max), abs(R44min)]);

[intB,~,~]=CalcCSRkick(beamline, -240, 0, 3e-5);


r56=r(5,6);
r26=abs(rmid(2,6));


if  r26<1e-4
    cdis=1;
else
    cdis=exp((abs(r26))/(1e-4));
end
f1=abs(intB(1));
f2=abs(intB(2));


Rmax=max(r(:));
if Rmax<4 
    c1 = 1;
else
    c1 = (abs(Rmax));
end



y = [abs(cdis*c1*f1)+cdis+c1-2, abs(cdis*c1*f2)+cdis+c1-2];

end