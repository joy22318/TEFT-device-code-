clc; clear; close all;

EG1= importdata('gate_voltage_current_10.3nm.mat');
EG2= importdata('gate_voltage_current_8.9nm.mat');
EG3= importdata('gate_voltage_current_7.52nm.mat');
EG4= importdata('gate_voltage_current_6.1nm.mat');

Vds_all = 0.05*ones(1,20);
Vgs_all = EG1.Vgs;
Ids_all = EG4.I;


selected_Vds = [0.05];   % you can add more: [0.05 0.1 0.2]


figure;
hold on; grid on; box on;

colors = lines(length(selected_Vds));
markers = {'o','s','d','^','v','>','<','p','h','x','+'};

for i = 1:length(selected_Vds)
    Vds_val = selected_Vds(i);

    idx = abs(Vds_all - Vds_val) < 1e-6;

    Vgs = Vgs_all(idx);
    Ids = Ids_all(idx);

    Ids(Ids <= 0) = 1e-18; % avoid log(0)

    semilogy(Vgs, Ids, ...
        'LineWidth', 1.8, ...
        'Color', colors(i,:), ...
        'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerSize', 6);
end

set(gca, 'YScale', 'log');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (A)');
title('I_{DS} vs V_{GS} (Log Scale)');
xlim([0 2]);

legendStrings = arrayfun(@(x) sprintf('V_{DS} = %.3f V', x), selected_Vds, 'UniformOutput', false);
legend(legendStrings, 'Location', 'best');

hold off;


fprintf('\n================ Subthreshold Swing (SS) Results ================\n');

for i = 1:length(selected_Vds)

    Vds_val = selected_Vds(i);
    idx = abs(Vds_all - Vds_val) < 1e-6;

    Vgs = Vgs_all(idx);
    Ids = Ids_all(idx);

    Ids(Ids <= 0) = 1e-16;

    logIds = log10(Ids);


    sub_idx = (Vgs > 0.2 & Vgs < 1.5); % set off current range

    Vgs_sub = Vgs(sub_idx);
    logIds_sub = logIds(sub_idx);

    % Sort data
    [Vgs_sub, sortIdx] = sort(Vgs_sub);
    logIds_sub = logIds_sub(sortIdx);

    coeffs = polyfit(Vgs_sub, logIds_sub, 1);
    slope = coeffs(1);  % d(logId)/d(Vgs)


    SS = -(1 / slope) * 1000;

    fprintf('Vds = %.3f V :  SS = %.2f mV/dec\n', Vds_val, SS);

    
    figure(100+i);
    semilogy(Vgs, Ids, 'b.-', 'LineWidth', 1.2);
    hold on; grid on; box on;

    Vgs_fit = linspace(min(Vgs_sub), max(Vgs_sub), 100);
    logIds_fit = polyval(coeffs, Vgs_fit);

    plot(Vgs_fit, 10.^logIds_fit, 'r--', 'LineWidth', 1.6);

    title(sprintf('Subthreshold Fit for V_{DS} = %.2f V', Vds_val));
    
    xlabel('V_{GS} (V)');
    ylabel('I_{DS} (A)');
   
    legend('Data', 'Linear Fit (Subthreshold)');
end


