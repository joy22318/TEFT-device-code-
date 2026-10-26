clc; clear; close all;

EG1 = importdata('gate_voltage_current_6.1nm.mat');
EG2 = importdata('gate_voltage_current_7.52nm.mat');
EG3 = importdata('gate_voltage_current_8.9nm.mat');
EG4 = importdata('gate_voltage_current_10.3nm.mat');

Vds_all = 0.05*ones(1,20);
Vgs_all = 1.26706*EG2.Vgs;

% --- FIX 1: Extract absolute raw currents (DO NOT apply log10 here) ---
% semilogy will handle the log scaling automatically for the plot.
Ids1 = abs(real(EG1.I'));
Ids2 = abs(real(EG2.I'));
Ids3 = abs(real(EG3.I'));
Ids4 = abs(real(EG4.I'));

Length = zeros(1,4);
Length(1) = EG1.Length_Device;
Length(2) = EG2.Length_Device;
Length(3) = EG3.Length_Device;
Length(4) = EG4.Length_Device;

figure('Color', [1 1 1]); 
hold on;
grid on;
box on;

% --- Plotting with semilogy ---
% Using custom markers and distinct colors to mimic the reference style
p1 = semilogy(Vgs_all, Ids1, 'b-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'); % Blue Square
p2 = semilogy(Vgs_all, Ids2, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'); % Red/Orange Circle
p3 = semilogy(Vgs_all, Ids3, 'y-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'); % Yellow Triangle
p4 = semilogy(Vgs_all, Ids4, 'm-v', 'LineWidth', 1.5, 'MarkerFaceColor', 'w'); % Purple/Magenta Triangle

% --- Generate legend labels dynamically ---
legend_labels = arrayfun(@(Length) sprintf('Length = %.2f nm', Length*0.1), Length, 'UniformOutput', false);

% --- Axis Labeling & Styling ---
% Match the X-axis label style from your reference image
xlabel('V_{TG} - V_{BG} (V)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('I_{DS} (A)', 'FontSize', 20, 'FontWeight', 'bold');

% --- Legend & Grid Enhancement ---
legend(legend_labels, 'Location', 'best', 'FontSize', 16); 
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5);
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on'); 
set(gcf, 'Position', [100 100 1000 700]); % Slightly taller ratio to match reference image
set(gca, 'FontSize', 18);
xlim([min(Vgs_all) max(Vgs_all)]);