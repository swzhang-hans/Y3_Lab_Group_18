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
