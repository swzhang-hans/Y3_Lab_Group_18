% Parameters
Re = 4e5; % Target Reynolds number
rho_air = 1.225; % Density of air (kg/m^3)
mu_air = 1.81e-5; % Dynamic viscosity of air (Pa·s)
root_chord = 0.225; 
rho_manometer = 866; % Density of manometer fluid (kg/m^3)
g = 9.81; 

% Calculate wind speed (U)
U = (Re * mu_air) / (rho_air * root_chord);

% Calculate pressure difference (Delta P)
Delta_P = 0.5 * rho_air * U^2;

% Calculate manometer height (h)
h = Delta_P / (rho_manometer * g);


fprintf('Required wind speed: %.2f m/s\n', U);
fprintf('Manometer reading: %.2f mm\n', h * 1000); % Convert to mm
