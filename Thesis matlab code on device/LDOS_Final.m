clc
clear 
close all

%% Load Data
LDOS2 = importdata('zigzag.mat');     % Zigzag
LDOS22 = importdata('armchair.mat');  % Armchair

x = -1:128;

% Fit LDOS Data
[fitresult, ~] = createFit(x, LDOS2);
[fitresult1, ~] = createFit(x, LDOS22);

% Refined x-axis for smooth plotting
x = -1:0.001:128;
yy = fitresult(x);      % Zigzag
yy1 = fitresult1(x);    % Armchair

% Generate y-axis scale in nm
xx = linspace(0, 5.82, length(x));

%% Plotting
figure()
hold on

% Plot Armchair with soft red dashed line
plot(xx, yy1, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 2)

% Plot Zigzag with blue solid line
plot(xx, yy, 'b', 'LineWidth', 2)

% Axis Labels
xlabel('y (nm)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 12)
ylabel('LDOS (arb. units)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 12)

% Legend
legend({'armchair', 'zigzag'}, ...
       'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 12, ...
       'Location', 'north')

% Axis Properties
set(gca, 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 12)
pbaspect([2 1 1])
xlim([xx(1), xx(end)])

% Optional: add grid for better readability
grid on
box on


% Move Left Edge further left
annotation('textbox', [0.16 0.65 0.05 0.05], ...
    'String', 'Edge', 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontName', 'Times', 'FontWeight', 'bold');

% Bulk remains unchanged
annotation('textbox', [0.45 0.65 0.05 0.05], ...
    'String', 'Bulk', 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontName', 'Times', 'FontWeight', 'bold');

% Move Right Edge further right
annotation('textbox', [0.78 0.65 0.05 0.05], ...
    'String', 'Edge', 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'FontName', 'Times', 'FontWeight', 'bold');

% Double-Headed Arrows
annotation('doublearrow', [0.13 0.27], [0.6 0.6], 'LineWidth', 1.5); % Left Edge
annotation('doublearrow', [0.27 0.77], [0.6 0.6], 'LineWidth', 1.5); % Right Edge
annotation('doublearrow', [0.77 .91], [0.6 0.6], 'LineWidth', 1.5); % Bulk (middle)
