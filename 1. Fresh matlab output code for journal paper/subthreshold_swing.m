clc
clear
close all
format long;
set(0,'DefaultAxesFontName', 'Latex');
set(0,'DefaultAxesFontSize', 12);
% %% Load Data
EG2= importdata('Ids_Vds_10nm_length.mat');
EG3= importdata('Ids_Vg_8.56nm_length.mat');
EG4= importdata('Ids_Vg_6.82nm_length.mat');
EG5= importdata('Ids_Vg_5.072nm_length.mat');
Vtb=EG2.Vtb;

SS(1)=Subthreshold(EG5);
SS(2)=Subthreshold(EG4);
SS(3)=Subthreshold(EG3);
SS(4)=Subthreshold(EG2);
Length = [6.1, 7.52, 8.9, 10.31]; % Define the channel lengths corresponding to the data
Subthreshold_swing = SS; % Assign the calculated subthreshold swings to the variable for plotting
figure;
grid on
set(gca,'FontSize',20)
hold on;  % This command will keep the existing plot and allow adding new ones

plot(Length, Subthreshold_swing, '-o', 'LineWidth', 2, 'Color', 'b');
ylim([0 350]);


xlabel('Channel length (nm)', 'Interpreter','latex')
ylabel('Subthreshold Swing (mv/dec)', 'Interpreter','latex')
