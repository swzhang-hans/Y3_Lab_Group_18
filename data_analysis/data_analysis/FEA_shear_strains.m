% FEA results
% ====================================
clear;clc

% Array Inputs for loads and shear strains of FEA
FEA.real.loads = [0	7.390417941	10	22.17125382	30	36.9520897	40];
FEA.real.gamma_leadingEdge = [] ;
FEA.real.gamma_frontSpar = [0	3.8581E-06	1.72059E-06	1.15743E-05	5.3433E-06	1.92905E-05	6.8823E-06];
FEA.real.gamma_rearSpar = [0	6.33045E-06	5.21565E-06	1.89914E-05	1.56472E-05	3.16523E-05	2.08632E-05];
FEA.real.gamma_skinTop = [0	8.84768E-06	1.19186E-05	0.000026543	3.57564E-05	4.42385E-05	4.76746E-05];
FEA.real.gamma_skinBottom = [0	1.27128E-06	1.20026E-05	3.81383E-06	3.6008E-05	6.35628E-06	4.80087E-05];

FEA.simple.gamma_leadingEdge = [];
FEA.simple.gamma_frontSpar = [];
FEA.simple.gamma_rearSpar = [];
FEA.simple.gamma_skinTop = [];
FEA.simple.gamma_skinBottom = [];

save ("FEA_shear_strains.mat","FEA");