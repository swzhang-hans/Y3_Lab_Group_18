clc
clear


A1=8478*10^(-6); %m^2
A2=24628*10^(-6);
G=28*10^(9); %Pa
t18e=0.64*10^(-3); %m
t18i=1.7*10^(-3);
t12=2.34*10^(-3);
t23=0.64*10^(-3);
t36=1.59*10^(-3);
t78=2.34*10^(-3);
t67=t23;
b18e=318*10^(-3);
b18i=81*10^(-3);
b12=17.4*10^(-3);
b23=286.66*10^(-3);
b36=81*10^(-3);
b78=b12;
b67=t23;
P=50; %lbs
P=P*0.453592*9.81; %N
L=b12+b23; %m
T=P*L; %T=P*L
syms q1 q2 theta_z
equs = [2*q1*A1+2*q2*A2==T, theta_z==1/(2*A1*G)*(q1*b18e/t18e+(q1-q2)*b18i/t18i) , theta_z==1/(2*A2*G)*(q2*(b23/t23+b36/t36)+(q2-q1)*b18i/t18i)];
vars = [q1,q2,theta_z];
[solq1,solq2,soltheta_z]=solve(equs,vars);

leading_edge=double(solq1/t18e*G);
front_spar=double((solq2-solq1)/t18i*G);
top_skin=double(solq2/t23*G);
rear_spar=double(solq2/t36*G);


%%
I_xx = 8.971428*10^-7;
E = 73.1*10^9;
L = 2;

Tip_deflection = (50 * L^3) / (3 * (I_xx) * (E));

sigma_max_root = (50 * 2000 * 40*10^-3) / I_xx;

L_eff = 4 * L;

P_crit = (pi ^ 2 * E * I_xx) / (L_eff)^2;





