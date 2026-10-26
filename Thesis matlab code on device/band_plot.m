clc 
close all
clear all

% Plot Wannier90 data (sky blue dots)
EF=3.13195
data=importdata('wte2_band.dat');
nbnd = 44;
kpoints = 201;
for i = 1:kpoints
    kwan(i) = data(i,1);
end

x = 1;
for j = 1:kpoints:kpoints*nbnd
    for i = 1:kpoints
        Ewan(x,i) = data(j+i-1,2) - EF;
    end
    x = x + 1;
end

p2 = plot(kwan, Ewan,'k', 'LineWidth', 1.5); % Store plot handle for red line

% Draw vertical lines (arrows without heads)
line([0.89813, 0.89813], [-15, 8.39784], 'Color', 'k');
line([1.79625, 1.79625], [-15, 8.39784], 'Color', 'k');
line([2.29191, 2.29191], [-15, 8.39784], 'Color', 'k');
line([0, 2.78756], [0, 0], 'Color', 'k');

% Set custom x-ticks
set(gca, 'XTick', [0.00000, 0.89813,1.79625]);
set(gca, 'XTickLabel', {'X', '$\Gamma$', 'X'}, 'TickLabelInterpreter', 'latex');
axis tight

% Set the y-axis and x-axis labels
ylabel('E-EF (eV)');
xlabel('Wave Vector');

% Set axis limits
xlim([0, 1.79625]);
ylim([-1.5, 1.5]);


hold off;

Bandgap=min(Ewan(29,:))-max(Ewan(28,:))