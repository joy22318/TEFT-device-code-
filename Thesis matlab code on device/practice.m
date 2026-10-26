clc;
close all;
clear all;

% Constants
Ef = 0;               % Fermi energy (eV)
Vds = 0.05;           % Voltage (V)
kT = 0.025;          % Thermal energy at room temperature (0.0259 eV at 300 K)

% Energy range
mu1 = Ef+Vds/2;
mu2 = Ef - Vds/2;
E = linspace(Ef-4*kT, Ef+4*kT, 100);  % Energy from mu2 to mu1
dE = E(3) - E(2);            % Energy step
delphi = 1;
Ite = 0;

% Fermi-Dirac distribution (Boltzmann approximation)
f1 = 1 ./ (1 + exp((E) / kT));
f2 = 1 ./ (1 + exp((E) / kT));

% Derivatives of f1 and f2
df1 = -f1 .* (1 - f1) / kT;
df2 = -f2 .* (1 - f2) / kT;

F_T = -(1/(4*kT)).*(sech((E-Ef)./(2*kT))).^2;

% Plotting
figure;
subplot(2,1,1);
plot(E, f1, '-o', 'DisplayName', 'f1 (Ef)');
hold on;
plot(E, f2, '-x', 'DisplayName', 'f2 (Ef - Vds)');
xlabel('Energy (eV)');
ylabel('f(E)');
title('Boltzmann Approximation vs Energy');
legend;
grid on;
hold off;

% Plotting the derivatives
subplot(2,1,2);
plot(E, F_T, '-o', 'DisplayName', 'df1/dE');
hold on;
plot(E, F_T, '-x', 'DisplayName', 'df2/dE');
xlabel('Energy (eV)');
ylabel('df/dE');
title('Derivative of Boltzmann Approximation');
legend;
grid on;
hold off;








% Define the range for x
x = linspace(-2*pi, 2*pi, 1000);

% Compute y = sec(x)
y_sec = sech(x);

% Compute y = sech(x)
y_sech = 2./(exp(x)+exp(-x));

% Plot sec(x) and sech(x)
figure;
hold on;

% Plot sec(x)
plot(x, y_sec, 'r', 'LineWidth', 2, 'DisplayName', 'sec(x)');

% Plot sech(x)
plot(x, y_sech, '--b', 'LineWidth', 2, 'DisplayName', 'sech(x)');

% Labeling the graph
title('Plot of (sech(x)) 2/(e^(x)+e^(-x))');
xlabel('x');
ylabel('y');
legend;
grid on;
hold off;
