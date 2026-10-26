clc
clear
close all
caption = 'wte2';
% I_vs_Vgs_temp_fixed.m  --  1T'-WTe2 TIFET, Ids-Vgs and Cq at several temperatures
%
% Fixes relative to I_vs_Vgs_temp.m
%  1) Energy grid now covers the full Fermi window [mu2 - 10kT, mu1 + 10kT].
%     The old grid E = linspace(mu2,mu1,10) cut the thermal tails, so the
%     300 K current was underestimated (~47 %) and the ON current looked
%     "doubled" at 4 K. The same truncation/under-sampling corrupted Cq(T).
%  2) Prefactor q^2/h (NOT 2q^2/h): the 44-orbital basis
%     (4 Te x 6 + 2 W x 10) is a spinor basis, so T(E) already counts spin.
%  3) kB = 8.617333262e-5 eV/K (old value 8.3333e-5 was 3.3 % low).
%  4) Temperature defined once (old script set Temp twice).
%  5) Coherent NEGF: T(E) and DOS(E) do not depend on temperature, so they
%     are computed once and every temperature is post-processed from them.
%     Fermi-window integrals are done analytically per energy bin, which is
%     exact for piecewise-constant T(E) and stays accurate at 4 K.
%  6) Unused Gn/Gp removed; transmission from corner block G_1N only;
%     Hamiltonian and lead self-energies computed once (not per Vgs).

%% ======================= USER SETTINGS =======================
Temps = [300 77 4];            % K, all temperatures handled in one run
Vds   = 0.05;                  % V
mu2   = 3.0;                   % drain electrochemical potential (eV)
mu1   = mu2 + Vds;             % source electrochemical potential (eV) = 3.05
Ef    = 3.015;                 % equilibrium Fermi level used for Cq (eV)
Vgs   = linspace(0,1.5,20);    % gate voltage (V)
dE    = 2e-3;                  % energy step (eV): must resolve T(E), not kT
nkT   = 10;                    % Fermi-tail margin in units of kT at max(Temps)

%% ======================= CONSTANTS (CODATA) =======================
q        = 1.602176634e-19;    % C
h_planck = 6.62607015e-34;     % J s
kb       = 8.617333262e-5;     % eV/K
I0       = q^2/h_planck;       % A/V, spin already inside T(E): no factor 2

%% ======================= GEOMETRY (unchanged) =======================
ZBNR = 12;                     % zigzag nanoribbon index: 4,8,12,16,...
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn;            % 3,6,9,12,...
Unit_Cells = 30;

xa(1)=0; xa(2)=0; xa(3)=0;
ya(1)=0; ya(2)=1; ya(3)=2;
for j = 0:num_cell-4
    xa(4+j) = xa(1+j)+0;
    ya(4+j) = ya(1+j)+3;
end
za = zeros(1,length(xa))';
main_unit = [xa', ya', za];
cp = [xa', ya', za];

xb(1)=1; xb(2)=1; xb(3)=1;
yb(1)=0; yb(2)=1; yb(3)=2;
for j = 0:num_cell-4
    xb(4+j) = xb(1+j)+0;
    yb(4+j) = yb(1+j)+3;
end
zb = zeros(1,length(xb))';
cpp = [xb', yb', zb];

Device_unit = [];
for i = 1:Unit_Cells
    for j = 1:size(main_unit,1)
        Device_unit(end+1,:) = main_unit(j,:) + [(i-1) 0 0];
    end
end

atom_factors = [0.00000044, 1.31473454, 6.35624338;
                1.74897061, 4.43335644, 5.73543738;
               -0.00000002, 5.44760823, 9.92929379;
                1.74897108, 2.22834683, 9.30776171;
                1.74897050, 6.06754902, 7.93340128;
                0.00000040, 3.81415493, 7.73232446];
unit_cell_length = [3.49794 6.33827 15.4319];

Device_atomic_position = [];
for i = 1:size(Device_unit,1)
    for j = 1:size(atom_factors,1)
        Device_atomic_position(end+1,:) = Device_unit(i,:).*unit_cell_length + atom_factors(j,:);
    end
end
x_uniq = unique(Device_atomic_position(:,1));
Length_Device = max(x_uniq) - min(x_uniq);     % Angstrom

na  = num_cell/3;
at1 = [1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6];
at  = at1;
for j = 1:na-1
    at = [at; at1];
end

%% ======================= HAMILTONIAN (built once) =======================
cell_pos = cp;  cell_pos1 = cpp;  atomst = at;  atomst1 = at;
proj = [6 6 6 6 10 10];
HNN  = 5;
Vext = 0;
eta  = 1e-6;
count1 = 44*num_cell;               % orbitals per slice
count2 = count1*Unit_Cells;         % orbitals in device

[bigH,bigU,bigUbar] = cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);
H = full(hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells));   % independent of E and Vgs

% Gate on-site potential (diagonal), same model as before:
% +lamda_v on Te1,Te2 (orbitals 1-12), -lamda_v on Te3,Te4 (13-24), 0 on W
NV = numel(Vgs);
vgate = zeros(count2,NV);
for kk = 1:NV
    lamda_v = Vgs(kk);
    d44 = zeros(44,1);
    d44(1:12)  = +lamda_v;
    d44(13:24) = -lamda_v;
    vgate(:,kk) = repmat(d44, num_cell*Unit_Cells, 1);
end

%% ======================= ENERGY GRID =======================
kTmax = kb*max(Temps);
E  = (mu2 - nkT*kTmax) : dE : (mu1 + nkT*kTmax);
NE = numel(E);
fprintf('Energy grid: %d points, %.4f to %.4f eV\n', NE, E(1), E(end));

%% ======================= NEGF: T(E,Vgs), DOS(E,Vgs) =======================
Tr  = zeros(NV,NE);                 % transmission (spin included)
DOS = zeros(NV,NE);                 % states/(eV*Angstrom)
B1  = 1:count1;                     % left-contact slice
BN  = count2-count1+1:count2;       % right-contact slice
dg  = 1:count2+1:count2^2;          % linear indices of the diagonal

for k1 = 1:NE
    [G1,~] = lead_eps(bigH,bigUbar,eta,E(k1),count1,count2);
    [G2,~] = lead_eps(bigH,bigU,eta,E(k1),count1,count2);
    sig1 = bigUbar*G1*bigU;         % left self-energy  (count1 x count1)
    sig2 = bigU*G2*bigUbar;         % right self-energy (count1 x count1)
    gam1 = 1i*(sig1 - sig1');
    gam2 = 1i*(sig2 - sig2');

    for kk = 1:NV
        A = -H;
        A(dg) = A(dg) + (E(k1) - vgate(:,kk)).';
        A(B1,B1) = A(B1,B1) - sig1;
        A(BN,BN) = A(BN,BN) - sig2;
        Gr  = inv(A);
        G1N = Gr(B1,BN);
        Tr(kk,k1)  = real(trace(gam1*G1N*gam2*G1N'));
        DOS(kk,k1) = -imag(trace(Gr))/(pi*Length_Device);
    end
    fprintf('E point %d / %d done\n', k1, NE);
end

%% ======================= TEMPERATURE POST-PROCESSING =======================
% Exact integral of the Fermi functions over each bin [E-dE/2, E+dE/2]:
%   w(k)  = int_bin [f(E-mu1) - f(E-mu2)] dE     (current window, eV)
%   wF(k) = int_bin F_T(E-Ef) dE, F_T = -df/dE    (Cq thermal kernel)
softplus = @(y) max(y,0) + log1p(exp(-abs(y)));     % log(1+e^y), overflow-safe
Fint     = @(x,mu,kT) -kT*softplus(-(x-mu)/kT);     % antiderivative of f
fermi    = @(x,mu,kT) 1./(1 + exp((x-mu)/kT));
Elo = E - dE/2;
Ehi = E + dE/2;

NT = numel(Temps);
I  = zeros(NT,NV);
Do = zeros(NT,NV);
for it = 1:NT
    kT = kb*Temps(it);
    w  = (Fint(Ehi,mu1,kT) - Fint(Elo,mu1,kT)) - (Fint(Ehi,mu2,kT) - Fint(Elo,mu2,kT));
    wF = fermi(Elo,Ef,kT) - fermi(Ehi,Ef,kT);
    fprintf('T = %3g K: sum(w) = %.5f eV (target %.5f), sum(wF) = %.5f (target 1)\n', ...
            Temps(it), sum(w), mu1-mu2, sum(wF));
    I(it,:)  = I0*(Tr*w.').';       % A
    Do(it,:) = (DOS*wF.').';        % states/(eV*Angstrom)
end

eV    = 1.602176634e-19;
A_ang = 1e-10;
Cq1 = (q^2/(eV*A_ang^2)).*Do;       % F/m^2
Cq2 = (q^2/(eV*A_ang)).*Do;         % F/m

%% ======================= SAVE (same file names as before) =======================
for it = 1:NT
    S1 = struct('Vgs',Vgs, 'I',I(it,:), 'Length_Device',Length_Device);
    save(sprintf('gate_voltage_current_temp%g.mat',Temps(it)), '-struct', 'S1');
    S2 = struct('Vgs',Vgs, 'Cq1',Cq1(it,:), 'Cq2',Cq2(it,:));
    save(sprintf('quantum_capacitance_data_temp%g.mat',Temps(it)), '-struct', 'S2');
end
save('TE_DOS_vs_E_Vgs.mat','E','Vgs','Tr','DOS','mu1','mu2','Ef','dE','Length_Device');

fprintf('\nBallistic limit for T = 2: 2(q^2/h)Vds = %.4f uA\n', 1e6*2*I0*Vds);
for it = 1:NT
    fprintf('T = %3g K: max(Ids) = %.4f uA, min(Ids) = %.4e A\n', ...
            Temps(it), 1e6*max(abs(I(it,:))), min(abs(I(it,:))));
end

%% ======================= PLOTS =======================
c_blue  = [0.000 0.447 0.741];
c_verm  = [0.835 0.369 0.000];
c_green = [0.000 0.620 0.451];
c_black = [0.100 0.100 0.100];
cols = {c_blue, c_verm, c_green};

figure('Color',[1 1 1]);
hold on; box on; grid on;
p = gobjects(1,NT);
for it = 1:NT
    p(it) = plot(Vgs, abs(I(it,:)), '-o', 'Color', cols{mod(it-1,3)+1}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log', 'FontSize', 18);
xlabel('V_{gs} (V)');
ylabel('I_{ds} (A)');
legend(p, arrayfun(@(t) sprintf('T = %g K',t), Temps, 'UniformOutput', false), 'Location', 'best');

figure('Color',[1 1 1]);
hold on; box on; grid on;
p_T = plot(E, Tr(1,:), '-', 'Color', c_black, 'LineWidth', 2);
p_mu1 = xline(mu1, 'LineStyle', '--', 'Color', c_verm, 'LineWidth', 1.5);
p_mu2 = xline(mu2, 'LineStyle', '--', 'Color', c_blue, 'LineWidth', 1.5);
set(gca, 'FontSize', 18);
xlabel('E (eV)');
ylabel('T(E)');
title(sprintf('Transmission at V_{gs} = %.2f V', Vgs(1)));

figure('Color',[1 1 1]);
hold on; box on; grid on;
p_c = gobjects(1,NT);
for it = 1:NT
    p_c(it) = plot(Vgs, Cq2(it,:), '-o', 'Color', cols{mod(it-1,3)+1}, 'LineWidth', 2);
end
set(gca, 'FontSize', 18);
xlabel('V_{tb} (V)');
ylabel('C_q (F/m)');
legend(p_c, arrayfun(@(t) sprintf('T = %g K',t), Temps, 'UniformOutput', false), 'Location', 'best');