%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low Speed Flow Past a High Aspect Ratio Swept Wing LAB - 2023/2024  %%%
%%% Academic supervisor: Dr. Oliver R. H. Buxton                        %%%
%%% Data processing script - based on the originally developed file by  %%%
%%% Emilia Juda, Aeronautics Class of 2018                              %%%
%%% Important notice- remember to close "data.xlsx" before running this %%%
%%% Last modification - Francisco Oliveira 08/11/2023                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc; close all;
%%% Some flags to run the code                                          %%%
NoRuns = 3;
NoSections = 3;
%%% For data interpolation purposes                                     %%%
N = 202; Nm=102;
RR = [1, Nm, N; N+1, Nm+N, 2*N; 1+2*N, Nm+2*N, 3*N];
%%% Describing the number of taps in each section:                      %%%
%%% Green: 15 | Yellow: 23 | Red: 23
Cond = [23, 23, 15];
%%% Data treatment: vectorial range per section                         %%%
Range = [1,Cond(1);Cond(1)+1,sum(Cond(1:2));sum(Cond(1:2))+1,sum(Cond)];
%%% Setting the location of x/c=0 tapping                               %%%
rmm = mean(Range,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Starting by dealing with the experimental results                   %%%
%%% Importing data from "data.xlsx"      
fprintf("\nImporting experimental data retrieved in the laboratory " + ...
    "session.\nRemember to close data.xlsx before running this script.\n")
raw_data = importdata('data.xlsx');
%%% Loading data acquired from the manometers for each section,         %%%
%%% and pitot tubes.                                                    %%%
raw_exp_data = raw_data.data(:, 3:end);
%%% aoa - tested angles of attack                                       %%%
aoa = raw_exp_data(73,:);
fprintf("\nTarget Reynolds number: %i.\n Achieved Reynolds " + ...
    "number: %.2f.\n\n", raw_exp_data(78,1),mean(raw_exp_data(81,:)));
%%% Importing tappings locations (as chord fractions)                   %%%
tap_loc = load('tapping_locations.txt');
tap_loc(any(isnan(tap_loc), 2), :) = [];
%%% Initializing Cp - matrix containing Cp for the tests run            %%%
Cp = zeros(sum(Cond), NoRuns);
%%% Extraction of pitot tube readings                                   %%%
% copy pitot values
pitot_1 = raw_exp_data(25:26, :);
pitot_2 = raw_exp_data(71:72, :);
%%% Computing static, dynamic and stagnation pressures - averaged       %%%
%%% between the two readings                                            %%%
p_static = 0.5*(pitot_1(1,:) + pitot_2(1,:));
p_stagn = 0.5*(pitot_1(2,:) + pitot_2(2,:));
p_dyn = 0.5*((pitot_1(2,:) - pitot_1(1,:)) + (pitot_2(2,:) - pitot_2(1,:)));
%%% Computing Cp with experimental data                                 %%%
A = raw_exp_data(1:68,:);
A(any(isnan(A), 2), :) = []; A(24:25,:)=[];
section_raw_data = raw_exp_data(1:68,:);
for j = 1:NoRuns 
    Cp(:,j) = (A(:,j)-p_stagn(j))./p_dyn(j);
end
%%% Setting chord lengths at different sections                         %%%
span = [115; 300; 531];
c = -0.2167 * span + 255;
c = c/1000; % converting to metres
%%% interpolate readings taken during the lab using pchip               %%%
xx = 0:0.01:1;
for i = 1:NoRuns
    for j = 1:NoSections 
        Cp_spline(i, RR(j,1):RR(j,2)-1) = pchip(tap_loc(Range(j,1):rmm(j)),Cp(Range(j,1):rmm(j), i), xx);
        Cp_spline(i, RR(j,2):RR(j,3)) = pchip(tap_loc(rmm(j):Range(j,2)),Cp(rmm(j):Range(j,2), i), xx);
        Cp_TE(i,j) = 0.5 * (Cp_spline(i, RR(j,2)-1) + Cp_spline(i, RR(j,3)));
    end
end
Cp_TE = transpose(Cp_TE);
%%% Updating Range to add the TE data                                   %%%
Range_new = [Range(1,1), Range(1,2)+2; Range(2,1)+2, Range(2,2)+4; Range(3,1)+4, Range(3,2)+6];
rmm = mean(Range_new, 2);
%%% Add the TE data to Cp matrix                                        %%%
tap_loc = [1; tap_loc(Range(1,1):Range(1,2)); 1; 1; tap_loc(Range(2,1):Range(2,2));...
    1; 1; tap_loc(Range(3,1):Range(3,2)); 1]; 
Cp = [Cp_TE(1, :); Cp(Range(1,1):Range(1,2), :); Cp_TE(1, :); Cp_TE(2, :); ...
    Cp(Range(2,1):Range(2,2), :); Cp_TE(2, :); Cp_TE(3, :); Cp(Range(3,1)...
    :Range(3,2), :); Cp_TE(3, :)];
%%% interpolate once again, this time with averaged trailing edge       %%%
for i = 1:NoRuns
    for j = 1:NoSections
        Cp_spline(i, RR(j,1):RR(j,2)-1) = pchip(tap_loc(Range_new(j,1):rmm(j)),Cp(Range_new(j,1):rmm(j), i), xx);
        Cp_spline(i, RR(j,2):RR(j,3)) = pchip(tap_loc(rmm(j):Range_new(j,2)),Cp(rmm(j):Range_new(j,2), i), xx);
    end
end
%%% Finished gathering Experimental data                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Importing panel code outputs                                        %%%
section_length = importdata('output_data/deg0/sect_001'); 
section_length = size(section_length(:,1),1);
rbpc = [1; 1 + (NoSections-2)*section_length; 1 + (NoSections-1)*section_length ];
repc = [section_length; section_length*(NoSections-1); section_length*NoSections];
%%% Cp_panel: [NoRuns, grouped sections data size, airfoil 2 sides]     %%%
Cp_panel = zeros(NoRuns, section_length*NoSections, 2);
%%% Load output from panel code                                         %%%
for i = 1:NoRuns
    str =['./output_data/deg' num2str(aoa(i))];
    for j = 1:NoSections
        filename = sprintf('%s/sect_00%s', str, num2str(j));
        Cp_panel(i, rbpc(j):repc(j), :) = importdata(filename);
    end
end
%%% Following the same sign convention as used in the experimental data %%%
Cp_panel(:, :, 2) = -Cp_panel(:, :, 2);
%%% Finished gathering Panel code output data                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computing sectional Lift coefficient Cl from panel code and         %%%
%%% experimental data output                                            %%%
for i = 1:NoRuns
    for j = 1:NoSections
        Cl(i,j) = - trapz(xx, Cp_spline(i, RR(j,1):RR(j,2)-1)) + trapz(xx, Cp_spline(i, RR(j,2):RR(j,3)));
        Cl_panel(i,j) = trapz(Cp_panel(i, rbpc(j):repc(j), 1), Cp_panel(i, rbpc(j):repc(j), 2));
    end
end
%%% Compute theoretical predition for Cl                                %%%
%%% See Fig. 19, reference [3] - handout                                %%%
%%% at 0 deg, Cl = 0.4667; at 4 deg, Cl = 0.9                           %%%
xx_t = 0:aoa(3);
Cl_theory = polyval([0.12 0.52], xx_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Starting to plot                                                    %%%
%%% Setting axes for Cl vs AoA plotting                                 %%%
cl_max = max(Cl_theory);
ax_cl = [0 aoa(3) 0 cl_max];
%%% Setting plotting style for Cp vs x/c graphs
crosses = {'xr', 'xy', 'xg'};
lines = {'-r', '-y', '-g'};
for fig = 1:NoRuns
    figure(fig)
    for j = 1:NoSections
        subplot(1,3,j)
        plot(tap_loc(Range(j,1):Range(j,2)), Cp(Range(j,1):Range(j,2), fig), crosses{j});
        hold all
        plot(xx, Cp_spline(fig, RR(j,1):RR(j,2)-1), lines{j});
        plot(xx, Cp_spline(fig, RR(j,2):RR(j,3)), lines{j});
        plot(Cp_panel(fig, rbpc(j):repc(j), 1), Cp_panel(fig, rbpc(j):repc(j), 2), '-k');
    end    
end
%%% Setting plotting style for Cl vs AoA plot                           %%%
lines2 = {'-xr', '-xy', '-xg'};
%%% Plotting experimental results                                       %%%
fig = fig +1;
figure(fig)
hold on
plot(xx_t, Cl_theory, '-k');
for i = 1:NoRuns
   plot(aoa, Cl(:,i), lines2{i});
end
axis(ax_cl);
grid on;
title('Experimental C_L-\alpha variations of the three sections')
legend('infinite wing','root section', 'middle section', 'tip section','Location','southeast');
xlabel('$\mathrm{angle\:of\:attack\:(degrees)}$','Interpreter','latex');
ylabel('$\mathrm{lift\:coefficient}\:\:C_l$','Interpreter','latex');
saveas(gcf,'exp_cl','epsc')

%%% Plotting panel code results                                         %%%
fig = fig +1;
figure(fig)
hold on
plot(xx_t, Cl_theory, '-k');
for i = 1:NoRuns
    plot(aoa, Cl_panel(:,i), lines2{i});
end
axis(ax_cl);
grid on;
title('Panel code C_L-\alpha predictions of the three sections')
legend('infinite wing', 'root section', 'middle section', 'tip section','Location','southeast');
xlabel('$\mathrm{angle\:of\:attack\:(degrees)}$','Interpreter','latex');
ylabel('$\mathrm{lift\:coefficient}\:\:C_l$','Interpreter','latex');
saveas(gcf,'panel_Cl','epsc')

%%% Compute the relative error in Cl prediction                         %%%
disp('Below are listed relative Cl errors. Every row corresponds to a different run, every column corresponds to a different section.');
Cl_error = abs((Cl - Cl_panel)./(Cl))
% row - run number; column - section number

disp('Select the best prediction');
prompt = 'Input run number (row number): ';
run_best = input(prompt);
prompt = 'Input section number (column number) [1 - red; 2 - yellow; 3 - green]: ';
section_best = input(prompt);

disp('Select the worst prediction');
prompt = 'Input run number (row number): ';
run_worst = input(prompt);
prompt = 'Input section number (column number) [1 - red; 2 - yellow; 3 - green]: ';
section_worst = input(prompt);

ax_pred = [0, 1, -inf, inf];

%%% Plotting best prediction                                            %%%
fig = fig + 1;
figure(fig)
h1 = plot(tap_loc(Range(section_best,1):Range(section_best,2)), -Cp(Range(section_best,1):Range(section_best,2), run_best), 'xk');
hold all
h2 = plot(xx, -Cp_spline(run_best, RR(section_best,1):(RR(section_best,2)-1)), '-k');
h3 = plot(xx, -Cp_spline(run_best, RR(section_best,2):RR(section_best,3)), '-k');
h4 = plot(Cp_panel(run_best, rbpc(section_best):repc(section_best), 1), -Cp_panel(run_best, rbpc(section_best):repc(section_best), 2), '--b');
legend([h1 h2 h4], {'Experimental datapoints', 'Data fit', 'Panel code prediction'}, 'Location', 'Southeast');
xlabel('$\mathrm{position\:along\:the\:chord}$','Interpreter','latex');
ylabel('$\mathrm{pressure\:coefficient}\:\:-C_p$','Interpreter','latex');
title('Best panel code C_p prediction juxtaposed with experimental values');
axis ij;
saveas(gcf,'best_pred','epsc')

%%% Plotting worst prediction                                            %%%
fig = fig + 1;
figure(fig)
h1 = plot(tap_loc(Range(section_worst,1):Range(section_worst,2)), -Cp(Range(section_worst,1):Range(section_worst,2), run_worst), 'xk');
hold all
h2 = plot(xx, -Cp_spline(run_worst, RR(section_worst,1):(RR(section_worst,2)-1)), '-k');
h3 = plot(xx, -Cp_spline(run_worst, RR(section_worst,2):RR(section_worst,3)), '-k');
h4 = plot(Cp_panel(run_worst, rbpc(section_worst):repc(section_worst), 1), -Cp_panel(run_worst, rbpc(section_worst):repc(section_worst), 2), '--b');
legend([h1 h2 h4], {'Experimental datapoints', 'Data fit', 'Panel code prediction'}, 'Location', 'Southeast');
xlabel('$\mathrm{position\:along\:the\:chord}$','Interpreter','latex');
ylabel('$\mathrm{pressure\:coefficient}\:\:-C_p$','Interpreter','latex');
title('Worst panel code C_p prediction juxtaposed with experimental values');
axis ij;
saveas(gcf,'worst_pred','epsc')
