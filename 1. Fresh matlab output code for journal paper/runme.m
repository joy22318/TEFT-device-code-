clc; clear; close all;

% EG1= importdata('gate_voltage_current_6.1nm.mat');
% EG2= importdata('gate_voltage_current_7.52nm.mat');
% EG3= importdata('gate_voltage_current_8.9nm.mat');
EG1= importdata('gate_voltage_current_edge_Vacancy_0.85_percent.mat');
EG2= importdata('gate_voltage_current_edge_Vacancy_2.0_percent.mat');
EG3= importdata('gate_voltage_current_random_Vacancy_0.85_percent.mat');
EG4= importdata('gate_voltage_current_random_Vacancy_2.0_percent.mat');
EG5= importdata('gate_voltage_current_random_Vacancy_5.0_percent.mat');
EG6= importdata('gate_voltage_current_10.3nm.mat');

Vds_all = 0.05*ones(1,20);
Vgs_all = 1.26706*EG2.Vgs;
Ids1 = real(EG1.I');
Ids2 = real(EG2.I');
Ids3 = real(EG3.I');
Ids4 = real(EG4.I');
Ids5 = real(EG5.I');
Ids6 = real(EG6.I');
Id1(14:20)=0;Ids2(11:20)=0; Ids3(14:20)=0; Ids4(11:20)=0; Ids5(10:20)=0;



figure('Color', [1 1 1]); % Set background to white
hold on;
grid on;
box on;

% --- Plotting with Distinguishable Features ---
% Red Solid line
% % % % % EG1= importdata('gate_voltage_current_edge_Vacancy_0.85_percent');
% % % % % EG2= importdata('gate_voltage_current_edge_Vacancy_2.0_percent');
% % % % % EG3= importdata('gate_voltage_current_random_Vacancy_0.85_percent');
% % % % % EG4= importdata('gate_voltage_current_random_Vacancy_2.0_percent');
% % % % % EG5= importdata('gate_voltage_current_random_Vacancy_5.0_percent');
% % % % % EG6= importdata('gate_voltage_current_10.3nm.mat');
% Green Dash-Dot line
p1 = plot(Vgs_all, Ids1, 'Color', '#D95319', 'LineStyle', '-',  'LineWidth', 1.5, 'Marker', 'o',  'DisplayName', '0.85% edge vacancy');
p2 = plot(Vgs_all, Ids2, 'Color', '#0072BD', 'LineStyle', '--', 'LineWidth', 1.5, 'Marker', 's',  'DisplayName', '2.0% edge vacancy');
p3 = plot(Vgs_all, Ids3, 'Color', '#EDB120', 'LineStyle', ':',  'LineWidth', 1.5, 'Marker', '^',  'DisplayName', '0.85% random vacancy');
p4 = plot(Vgs_all, Ids4, 'Color', '#7E2F8E', 'LineStyle', '-.', 'LineWidth', 1.5, 'Marker', 'v',  'DisplayName', '2.0% random vacancy');
p5 = plot(Vgs_all, Ids5, 'Color', '#2CA02C', 'LineStyle', '-',  'LineWidth', 1.5, 'Marker', 'd',  'DisplayName', '5.0% random vacancy');
p6 = plot(Vgs_all, Ids6, 'Color', '#000000', 'LineStyle', '-',  'LineWidth', 1.5, 'Marker', 'none', 'DisplayName', 'No vacancy');

% --- Optional: Add markers if the data points are sparse ---
% set(p1, 'Marker', 'o', 'MarkerIndices', 1:10:length(ku1));
% set(p2, 'Marker', 's', 'MarkerIndices', 1:10:length(ku1));
% set(p3, 'Marker', '^', 'MarkerIndices', 1:10:length(ku1));

%%%%%%%% xlabel and ylabel
xlabel('$V_{TG} - V_{BG}\;(\mathrm{V})$','Interpreter','latex', 'FontSize', 20)
ylabel('Current Ids (A) →', 'FontSize', 20)
% % % title('Distinguishable Magnetization Plot', 'FontSize', 14);

% --- Legend & Grid Enhancement ---
legend('show', 'Location', 'best', 'FontSize', 20);
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5);
set(gcf, 'Position', [100 100 1000 400]);
set(gca,'FontSize',18)
xlim([min(Vgs_all) max(Vgs_all)]);