%basic house cleaning
clear
clc
close all

%%
%set input
runtime="      "           ;               %e.g."2024-02-16_16-46-38" (a list of number on name of auto-output csv,repersent time)
Atmospheric_pressure=       ;             %(bar)
Atmospheric_temperature=       ;          %(degrees)

%%
%use the function given *(Processing function – supersoniclab2024.m) to analyze the data
[x,h,t,p0,p,xt,ts,prf]=supersoniclab2024(runtime, Atmospheric_pressure, Atmospheric_temperature); 
% OUTPUT(as the function itself said):
%  x     - Nozzle geometry, streamwise direction from nozzle part LE  [mm]
%  h     - Nozzle geometry, height of channel - h = f(x)              [mm]
%  t     - Test time                                                  [s]
%  p0    - Settling chamber total pressure, corrected for local flow  [Pa]
%  p     - Pressure matrix (time series x tap location)               [Pa]
%  xt    - Pressure tap location along centerline                     [mm]
%  ts    - Time limits for profiles (profile no. x [min(t) max(t)])   [s]
%  prf   - Normalised mean pressure profiles, p/p0                    [Pa]

%%
%caculate inviscid pressure rises in supersonic test cases we would expect to see in theory
%

%%
%plot other needed graphes (Normalised pressure ratio time series plot(1)&Normalised pressure profile plots(2) are already given) 

%Using 'stream-wise variation of normalised static pressure through the tunnel for each value of stagnation pressure tested'(2? not quite sure here) plot estimate the normal shock strength
%Annotate your pressure graph with the expected inviscid pressure rises for the transonic (based on your previous estimate of shock strength)



%%
%calculation
gData = importdata('icswtgeometryMach2.csv');
geo=gData.data(:,:);
%show geometry
x=geo(:,1);
h=geo(:,2);
h_u=h/2;
h_l=-h/2;
figure
plot(x,h_u,LineWidth=2,Color="b")
hold on
plot(x,h_l,LineWidth=2,Color="b")
axis('equal')
hold off

%Calculate Mach number with num-method
for i=1:length(h)
 area_ratio = h(i)/min(h);      
 f = @(M) (1/M) * 0.5787 * (1 + 0.2 * M^2)^3 - area_ratio; 

 M_sub(i) = fzero(f, [0.0001, 0.9999]);

 M_super(i) = fzero(f, [1.0001, 10]);
end
M1=find(h==min(h));

%Calculate normalised pressure
PperP0_1=(1+0.4/2*M_sub.^2).^(-1.4/0.4);
PperP0_2=(1+0.4/2*M_super.^2).^(-1.4/0.4);
%repair 
Prfsub2sup=[PperP0_1(1:M1(1)),PperP0_2(M1(2):end)];
Prfsup2sub=[PperP0_2(1:M1(1)),PperP0_1(M1(2):end)];
Prfsub2sub_special=PperP0_1;
figure
%subsonic to supersonic case
plot(x,Prfsub2sup,LineWidth=2,Color='r')
hold on
%supersonic to subsonic case
plot(x,Prfsup2sub,LineWidth=2,Color='b')
%supersonic to subsonic case with M at A* equals 1
plot(x,Prfsub2sub_special,LineWidth=1,LineStyle=":",Color='y')
%extra case where in all stage the flow is subsonic
M1=0.25;
for i=1:length(h)
 area_ratio = h(i)/h(1);      
 f = @(M) (1/M) * 0.5787 * (1 + 0.2 * M^2)^3/((1/M1) * 0.5787 * (1 + 0.2 * M1^2)^3) - area_ratio; 

 M_Sub2(i) = fzero(f, [0.0001, 0.9999]);
end
Prfsub2sub=(1+0.4/2*M_Sub2.^2).^(-1.4/0.4);
plot(x,Prfsub2sub,LineWidth=1,Color='k')

legend('subsonic to supersonic case','supersonic to subsonic case', ...
    'supersonic to subsonic case with M at A* equals 1','supersonic to subsonic case normal')
