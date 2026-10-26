clc;
clear;
close all;
format long;
set(0, 'DefaultAxesFontName', 'Latex');
set(0, 'DefaultAxesFontSize', 12);

% Import the data files
EG1 = importdata('gate_voltage_current_10.3nm.mat');
EG2 = importdata('gate_voltage_current_7.52nm.mat');

% Extract current and voltage data
I(:,1) = EG1.I;
I(:,2) = EG2.I;
Vtb = 1.26706*EG1.Vgs; % Voltage data

% Store device lengths
Length(1) = EG1.Length_Device;
Length(2) = EG2.Length_Device;

% Calculate transconductance for each current dataset (using diff to approximate the derivative)
g_m1 = diff(I(:,1)) ./ diff(Vtb);  % Transconductance for first dataset (10.32 nm length)
g_m2 = diff(I(:,2)) ./ diff(Vtb);  % Transconductance for second dataset (5.072 nm length)

% Create a new figure for plotting all transconductance together
figure;
hold on;

% Plot all transconductance against the voltage (Vtb)
h1 = plot(Vtb(2:end), g_m1(:,1), '-s', 'LineWidth', 2, 'Color', [0, 0, 1], 'DisplayName', sprintf('Length = %.2f nm', Length(1)*0.1));
h2 = plot(Vtb(2:end), g_m2(:,2), '-<', 'LineWidth', 2, 'Color', [0.5, 0, 0.5], 'DisplayName', sprintf('Length = %.2f nm', Length(2)*0.1));

% Finish plotting
hold off;

% Add a legend with unique labels for the different lengths
legend([h1, h2], {'Length = 10.31 nm', 'Length = 7.52 nm'}, ...
    'Location', 'northwest', 'FontSize', 18);


% Label the axes
xlabel('$V_{TG} - V_{BG} \; (\mathrm{V})$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Transconductance $g_m$ (S)', 'Interpreter', 'latex', 'FontSize', 20);

ax = gca;
ax.XAxis.FontSize = 16; % Makes X-ticks larger
ax.YAxis.FontSize = 16; % Keeps Y-ticks smaller

xlim([min(Vtb) max(Vtb)]);