clc
clear
close all
format long;

% Set up LaTeX plotting defaults
set(0,'DefaultAxesFontName', 'Latex');
set(0,'DefaultAxesFontSize', 12);

%% Load Data
% load('quantum_capacitance_data_temp300.mat')
% load('quantum_capacitance_data_temp77.mat')
% load('quantum_capacitance_data_temp4.mat')
load('quantum_capacitance_data_edge_Vacancy_0.85_temp300.mat')
% load('quantum_capacitance_data_edge_Vacancy_0.85_temp77.mat')
% load('quantum_capacitance_data_edge_Vacancy_0.85_temp4.mat')
% load('quantum_capacitance_data_edge_Vacancy_2.0_temp300.mat')
% load('quantum_capacitance_data_edge_Vacancy_2.0_temp77.mat')
% load('quantum_capacitance_data_edge_Vacancy_2.0_temp4.mat')

%%%%%%% High-dielectric insulator parameters %%%%%%%
eps_0 = 8.854e-12;          % Permittivity of free space (F/m)
eps_hfo2 = 25;              % Dielectric constant of HfO2
eps = eps_0 * eps_hfo2;     % Absolute permittivity of insulator

% Geometric dimensions (converted from nm to meters)
b = 5.503 * 1e-9;           % Device width in meters
d = 1 * 1e-9;               % Oxide thickness in meters
L = 10.31 * 1e-9;          % Channel length in meters
Width = b;                  % Alias variable for width (meters)

% Handle voltage variable naming compatibility
Vtb = 1.26706*Vgs;                  

% Layer data 
dTe2 = 5.73543738;
dTe3 = 9.92929379;
thfo2 = 10;                 % HfO2 physical thickness in nm
twte2 = (dTe3 - dTe2);      
eps_wte2 = 20;

%%%% Accurate EOT Calculation (Relative to SiO2 where eps_sio2 ~ 3.9)
eps_sio2 = 3.9;
EOT = thfo2 * (eps_sio2 / eps_hfo2); 

%%%% Capacitance Calculations (in SI units: Farads and Farads/meter)
% 1. Absolute Oxide Gate Capacitance (Farads)
C1 = (eps * L * b) / d;         

% 2. Oxide capacitance per unit width (F/m)
% For a double-gate (DG) structure, we multiply by 2 because gates act in parallel
Ci_per_width = 2 * (eps * L) / d; 

% Create an array of oxide capacitances matching the size of Vtb
Ci = Ci_per_width * ones(1, length(Vtb));  

% 3. Total series capacitance per unit width (F/m)
C_T = 1 ./ ((1 ./ Ci) + (1 ./ Cq2));    

%%%% Convert units from F/m to fF/um for clean plotting
% (F/m) * (1e15 fF / 1 F) * (1 m / 1e6 um) = (F/m) * 1e9 = fF/um
Cq_fFum = Cq2 * 1e9;     
C_T_fFum = C_T * 1e9;

%% Plotting Config
figure('Color','w');
plot(Vtb, Cq_fFum, '-or', ...
     Vtb, C_T_fFum, '-xk', ...
     'LineWidth', 1.3, 'MarkerSize', 6);

xlabel('$V_{TG} - V_{BG}$ (V)', 'Interpreter','latex', 'FontSize', 20);
ylabel('Capacitance (fF/$\mu$m)', 'Interpreter','latex', 'FontSize', 20);
legend({'Quantum Capacitance ($C_q$)', 'Total Capacitance ($C_T$)'}, ...
       'Interpreter','latex', 'FontSize', 16, 'Location', 'best');

set(gca, 'FontName', 'Times', 'FontSize', 16, 'LineWidth', 1);
grid off;
box on;
xlim([min(Vtb) max(Vtb)]);
% Save the figure
% print(gcf, 'f15.png', '-dpng', '-r600')  % 600 DPI PNG

%%%% Total average capacitance, Delay Time, and DED Calculations %%%
% Given Performance Parameters
Vds  = 0.05;                  % Drain-source / Supply voltage (V)
Vdd  = Vds;                   % Vdd definition matching standard formulas
Ion  = 3.67e-6;               % On-current in Amperes (7.72 uA -> 7.72e-6 A)
Ioff = 0.099e-12;            % Off-current in Amperes (pA range)
L_nm = 10.31e-9;                % Device length in nm

% 1. Compute Width in micrometers from the metric base
Width_um = Width * 1e6;       % meters * 1e6 = micrometers (~0.007503 um)

% 2. Calculate Normalized Ion (for engineering reference verification)
Ion_mum = Ion / Width_um;     % On-current normalized per micrometer width (A/um)

% 3. Extract Average Total Capacitance per unit width (fF/um)
C_avg_fFum = mean(C_T_fFum);  

% 4. Convert to absolute Capacitance Ctot (in Farads)
% (fF/um) * um = fF -> scaled by 1e-15 = absolute Farads
Ctot_absolute = C_avg_fFum * Width_um * 1e-15; 

%% Precise Metrics Calculation (SI Unit Standardized Verification)

% Delay Time (seconds) -> Convert to picoseconds (ps)
time_delay_sec = (Ctot_absolute * Vdd) / Ion;
time_delay_ps  = time_delay_sec * 1e12;

% Dynamic Energy Dissipation (Joules) -> Convert to attojoules (aJ)
% Reference Image Formula: DED = Vdd * Ion * tau_d
DED_joules = Vdd * Ion * time_delay_sec;
DED_aJ     = DED_joules * 1e18; 

% Dynamic Power Dissipation (Watts) over an active cycle
DPD_watts = 0.5 * C_avg_fFum * Vdd^2; 

% Static Power Dissipation (Watts)
SPD_watts = Vdd * Ioff;

% Product of DED and Time Delay (Joules * seconds)
DED_delay_product = DED_joules * time_delay_sec;

%% Displaying Results precisely in Command Window
fprintf('\n=================== RESULTS ===================\n');
fprintf('Device Width          : %.4f nm (%.6f um)\n', Width*1e9, Width_um);
fprintf('Avg Total Capacitance : %.4e F (%.4f fF/um)\n', Ctot_absolute, C_avg_fFum);
fprintf('Normalized Ion        : %.4f uA/um\n', Ion_mum * 1e6);
fprintf('-----------------------------------------------\n');
fprintf('Delay Time (tau_d)    : %.4f ps\n', time_delay_ps);
fprintf('Dynamic Energy (DED)  : %.4e aJ\n', DED_aJ);
fprintf('Dynamic Power (DPD)   : %.4e fj/um\n', DPD_watts);
fprintf('Static Power (SPD)    : %.4e W\n', SPD_watts);
fprintf('DED * Delay Product   : %.4e J*s\n', DED_delay_product);
fprintf('===============================================\n');