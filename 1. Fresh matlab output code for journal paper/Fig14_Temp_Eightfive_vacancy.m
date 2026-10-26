clc; clear; close all;

EG1 = importdata('gate_voltage_current_edge_Vacancy_0.85_percent.mat');
EG2 = importdata('gate_voltage_current_edge_Vacancy_0.85_temp77.mat');
EG3 = importdata('gate_voltage_current_edge_Vacancy_0.85_temp4.mat');


Vds_all = 0.05*ones(1,20);
Vgs_all = 1.26706*EG2.Vgs;

% --- FIX 1: Extract absolute raw currents (DO NOT apply log10 here) ---
% semilogy will handle the log scaling automatically for the plot.
Ids1 = abs(real(EG1.I'));
Ids2 = abs(real(EG2.I'));
Ids3 = abs(real(EG3.I'));


Temp = zeros(1,3);
Temp(1) = 300;
Temp(2) = 77;
Temp(3) = 4;


figure('Color', [1 1 1]); 
hold on;
grid on;
box on;

% --- 9 distinct colors ---
c1 = [0.0000 0.4470 0.7410];   % blue
c2 = [0.8500 0.3250 0.0980];   % vermillion
c3 = [0.4660 0.6740 0.1880];   % green
c4 = [0.4940 0.1840 0.5560];   % purple
c5 = [0.9290 0.6940 0.1250];   % orange
c6 = [0.0000 0.6000 0.6000];   % teal
c7 = [0.8000 0.2000 0.5500];   % magenta
c8 = [0.5500 0.3500 0.1500];   % brown
c9 = [0.1000 0.1000 0.1000];   % near-black
% --- Plotting with semilogy ---
% Using custom markers and distinct colors to mimic the reference style
p4 = semilogy(Vgs_all, Ids1, '-s', 'Color', c4, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerEdgeColor', c4, 'MarkerFaceColor', 'w'); % Purple Triangle-down
p5 = semilogy(Vgs_all, Ids2, '-o', 'Color', c5, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerEdgeColor', c5, 'MarkerFaceColor', 'w'); % Orange Diamond
p6 = semilogy(Vgs_all, Ids3, '-^', 'Color', c6, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerEdgeColor', c6, 'MarkerFaceColor', 'w'); % Teal Pentagram

% --- Generate legend labels dynamically ---
legend_labels = arrayfun(@(Temp) sprintf('T = %.0f K ', Temp), Temp, 'UniformOutput', false);

% --- Axis Labeling & Styling ---
% Match the X-axis label style from your reference image
xlabel('V_{TG} - V_{BG} (V)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('I_{DS} (A)', 'FontSize', 20, 'FontWeight', 'bold');
title('0.85 percent vacancy in channel')
% --- Legend & Grid Enhancement ---
legend(legend_labels, 'Location', 'best', 'FontSize', 16); 
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5);
set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'on'); 
set(gcf, 'Position', [100 100 1000 700]); % Slightly taller ratio to match reference image
set(gca, 'FontSize', 18);
xlim([min(Vgs_all) max(Vgs_all)]);