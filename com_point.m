function [y, cons] = com_point(x)



cons = [];                 % 约束初始化

beam_file='fvast2.dat';
EA=500;
E0=10;
L=0;
n_slices=100;
phi01=x(1);
phi02=x(2);
r1=x(3);
r2=x(4);
T566_1=x(5);
T566_2=x(6);
U5666_1=x(7);
U5666_2=0; x(8);



[energy_spread, rms_z, ~, kurtosis,minp]=longmat(beam_file, EA, E0, L, n_slices, phi01, phi02, r1, r2, T566_1, T566_2, U5666_1, U5666_2,0, 0);


if minp<-0.01

    c2=exp(abs(minp)*10)+5;

else

    c2=1;

end


if energy_spread<0.0015

    c1=1;

else

    c1=exp(2*energy_spread/0.0015);
end

f1=1e6*abs(rms_z-10e-6);

f2=kurtosis;

y=[c1*f1+c1+c2-2,c1*f2+c1+c2-2];

end