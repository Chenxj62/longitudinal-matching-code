function [intB,R56,beamline,l3]=schicanecalcu(r,L1,L2,L3,D1,D2,h,is_round)


%%l1->l2->l3

Lrange=[L1,L2,L3];
if is_round==1
Lrange=round(Lrange,2);
end

a1=2*asin(Lrange(1)/r/2);
a2=2*asin(Lrange(2)/r/2);
a3=2*asin(Lrange(3)/r/2);




l1=D1;
l2=D2;
l3=(1/2) * cot(a3) * (-4 * r + ...
    cos(a3) * (4 * r * cos(a2/2)^2 * sec(a1)^2 * sec(a2) - ...
    2 * sec(a1) * (2 * r + l1 * tan(a1)) - ...
    r * (1 + sec(a2))^2 * (-1 + tan(a1)^2) + ...
    2 * l2 * sec(a2) * tan(a2) + ...
    r * (-1 + tan(a1)^2) * tan(a2)^2));
R56=1/32*(-29*l1 - 18*l2 - 32*((2*(a1 + a2) + a3)*r + a3*r) + ...
    29*l1*sec(a1)^2 + 18*l2*sec(a2)^2 - 32*r*(-2 + cos(a3))*sin(a3) + ...
    32*r*cos(a3)*sin(a3) + 64*r*tan(a1) + 3*l1*tan(a1)^2 - ...
    32*sec(a1)*sin(a3)*(2*r + l1*tan(a1)) + 64*r*tan(a2) + ...
    14*l2*tan(a2)^2 + 32*sec(a2)*sin(a3)*(2*r + l2*tan(a2)));

ELE.B1 = createElement('dipole', a1*r, 1/r, 'Bchi1');
ELE.B2 = createElement('dipole',(a1+a2)*r, -1/r, 'Bchi2');
ELE.B3 = createElement('dipole', (a2+a3)*r, 1/r, 'Bchi3');
ELE.B4 = createElement('dipole', a3*r, -1/r, 'Bchi4');

E2_1=createElement('edge', 0, 1/r,a1, ' ');
E1_2=createElement('edge', 0, -1/r,-a1, ' ');
E2_2=createElement('edge', 0, -1/r,-a2, ' ');
E1_3=createElement('edge', 0, 1/r,a2, ' ');
E2_3=createElement('edge', 0, 1/r,a3, ' ');
E1_4=createElement('edge', 0, -1/r,-a3, ' ');

ELE.D1 = createElement('drift', l1, 0, 'Dchi1');
ELE.D2 = createElement('drift', l2, 0, 'Dchi2');
ELE.D3 = createElement('drift', l3, 0, 'Dchi3');

beamline=[ELE.B1,E2_1,ELE.D1,E1_2,ELE.B2,E2_2,ELE.D2,E1_3,ELE.B3,E2_3,ELE.D3,E1_4,ELE.B4];

beamline1=[ELE.B4,E1_4,ELE.D3,E2_3,ELE.B3,E1_3,ELE.D2,E2_2,ELE.B2,E1_2,ELE.D1,E2_1,ELE.B1];

[intB,~,~]=CalcCSRkick(beamline1,h,0,1);

end