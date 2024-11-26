clear 
clc
A1=8478;%mm^2
A2=24628;
G=28;%Gpa
t18e=0.64;%mm
t18i=1.7;
t12=2.34;
t23=0.64;
T36=1.59;
b18e=318;
b18i=81;
b12=17.4;
b23=286.66;
b36=81;
T=1;%T=P*L
syms q1 q2 theta_z
equs = [2*q1*A1+2*q2*A2==T, theta_z==1/(2*A1*G)*(q1*b18e/t18e+(q1-q2)*b18i/t18i) , theta_z==1/(2*A2*G)*(q2*(2*b23/t23+b36/t36)+(q2-q1)*b18i/t18i)];
vars = [q1,q2,theta_z];
[solq1,solq2,soltheta_z]=solve(equs,vars);

