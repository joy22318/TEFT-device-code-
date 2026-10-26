%% Vgs_Ids_edge_vacancy_T.m
% 1T'-WTe2 zigzag-nanoribbon TIFET: I_DS vs (V_TG - V_BG) with EDGE vacancies
% Temperatures: 300, 77, 4 K  |  edge vacancy: 0 (reference), 0.85, 2.0 %
% Coherent NEGF solved with the recursive Green's function (RGF) algorithm.
% Same Wannier Hamiltonian, leads and gate model as Vgs_Ids_vacancy.m.
% Needs on the path: cal_H_U_Ubar.m, hamiltonian_TMD_device_H.m, lead_eps.m
%
% Changes with respect to Vgs_Ids_vacancy.m
%  1) Vacancy -> orbital map fixed. The old map removed Te1 orbitals for a Te4
%     vacancy and Te1 + 4/6 of Te2 for a W1 vacancy; a deleted atom now loses
%     all of its own spin-orbitals (Kramers pairs kept intact).
%  2) Coherent transport => T(E) does not depend on temperature. T(E) is solved
%     once per vacancy configuration; temperature enters only through the
%     contact Fermi functions (and the sech^2 kernel of C_Q). The SAME vacancy
%     configuration is therefore used at 300, 77 and 4 K.
%  3) Energy grid spans [mu_D - n_kT*kB*Tmax, mu_S + n_kT*kB*Tmax]. The old
%     grid [mu_D, mu_S] kept only ~43 % of the 300 K Fermi window.
%  4) The 44-orbital SOC basis is spinful -> Landauer prefactor q^2/h (not 2q^2/h).
%  5) kB = 8.617333e-5 eV/K and kT = kB*T everywhere (old: 8.3333e-5 and a
%     fixed kT = 0.0259 eV in f1, f2 and F_T; Temp was also assigned twice).
%  6) Fixed RNG seed, nested configurations (0.85 % set is a subset of the
%     2.0 % set), no vacancies in the slices that touch the contacts.
%  7) Full 11880x11880 inversion and unused Gn, Gp replaced by RGF;
%     lead self-energies computed once per energy; trapz instead of sum.
clc; clear; close all;
caption = 'wte2';

%% ===================== User settings =====================
T_list      = [300 77 4];               % K
vac_list    = [0 0.85 2.0];             % edge vacancy, % of all atoms (0 = pristine reference)
vac_tag     = {'0','0.85','2.0'};       % used in the file names
seed        = 2026;                     % rerun with other seeds to average configurations
n_edge_rows = 2;                        % outermost atomic rows per edge that can lose atoms
n_buffer    = 2;                        % pristine slices kept next to each contact
U_vac       = 1e3;                      % eV, on-site energy of deleted orbitals

ZBNR        = 12;                       % must match the pristine run (12 -> 9 cells, W = 5.7 nm)
Unit_Cells  = 30;                       % slices along x (L = 10.3 nm)
lambda_list = linspace(0,1.5,20);       % staggered potential lambda_v (eV) = old "Vgs"

Vds  = 0.05;                            % V
muS  = 3.05;                            % eV, source (lead 1)
muD  = muS - Vds;                       % eV, drain  (lead 2) = 3.00
Ef   = 3.015;                           % eV, equilibrium Fermi level for C_Q (midpoint of window = 3.025)
dE   = 2.5e-3;                          % eV, Vds/dE must be an integer
n_kT = 10;                              % Fermi-tail margin in units of kB*T
eta  = 1e-6;                            % eV

t_ox = 14; eps_ox = 25; eps_wte2 = 1;   % HfO2 thickness (Angstrom) and permittivities

compute_DOS  = true;                    % false: ~40 % faster, no C_Q
validate_rgf = true;                    % one-shot check: RGF vs full-matrix solution
quick_test   = false;                   % true: coarse run to check the pipeline first
resume       = true;                    % continue from the checkpoint if it matches
ckpt_file    = 'edge_vacancy_checkpoint.mat';

if quick_test, lambda_list = linspace(0,1.5,4); dE = 1e-2; end

%% ===================== Constants =====================
q    = 1.602176634e-19;                 % C
h_pl = 6.62607015e-34;                  % J s
kB   = 8.617333262e-5;                  % eV/K
IE   = q^2/h_pl;                        % A/V; spin is explicit in the SOC basis

%% ===================== Geometry =====================
nn       = (ZBNR-4)/4 + 1;
num_cell = ZBNR - nn;                   % cells across the width
ax = 3.49794;  by = 6.33827;  cz = 15.4319;
unit_cell_length = [ax by cz];

main_unit     = [zeros(num_cell,1), (0:num_cell-1)', zeros(num_cell,1)];   % one slice
neighbor_unit = [ones(num_cell,1),  (0:num_cell-1)', zeros(num_cell,1)];   % next slice

atom_factors = [0.00000044, 1.31473454, 6.35624338;     % 1 Te  (bottom Te plane)
                1.74897061, 4.43335644, 5.73543738;     % 2 Te  (bottom Te plane)
               -0.00000002, 5.44760823, 9.92929379;     % 3 Te  (top Te plane)
                1.74897108, 2.22834683, 9.30776171;     % 4 Te  (top Te plane)
                1.74897050, 6.06754902, 7.93340128;     % 5 W
                0.00000040, 3.81415493, 7.73232446];    % 6 W
norb = [6 6 6 6 10 10];                 % spin-orbitals per atom (Te: p x 2, W: d x 2)
off  = [0 cumsum(norb(1:end-1))];       % orbital offset of each atom inside a cell
nob  = sum(norb);                       % 44 orbitals per cell

n_cell_dev  = num_cell*Unit_Cells;
Device_unit = zeros(n_cell_dev,3);
for i = 1:Unit_Cells
    Device_unit((i-1)*num_cell+(1:num_cell),:) = main_unit + [(i-1) 0 0];
end
n_at = 6*n_cell_dev;
Device_atomic_position = zeros(n_at,3);
for c = 1:n_cell_dev
    Device_atomic_position(6*(c-1)+(1:6),:) = Device_unit(c,:).*unit_cell_length + atom_factors;
end
Length_Device = max(Device_atomic_position(:,1)) - min(Device_atomic_position(:,1));  % Angstrom
L_dev = Unit_Cells*ax;                  % Angstrom, length of the N slices (DOS normalisation)
W_rib = num_cell*by;                    % Angstrom, ribbon width (C_Q per area)

t_wte2      = atom_factors(3,3) - atom_factors(2,3);        % 4.19386 A (outer Te-Te)
gate_factor = (2*t_ox*eps_wte2/eps_ox + t_wte2)/t_wte2;     % 1.26706 for t_ox = 14 A
fprintf('W = %.2f nm, L = %.2f nm, %d atoms, (V_TG-V_BG)/lambda = %.5f\n', ...
        W_rib/10, Length_Device/10, n_at, gate_factor);

%% ===================== Pristine Hamiltonian =====================
proj = [6 6 6 6 10 10];  HNN = 5;  Vext = 0;
atomst = repmat(1:6, num_cell, 1);
[bigH,bigU,bigUbar] = cal_H_U_Ubar(caption,main_unit,neighbor_unit,atomst,atomst,proj,HNN,Vext);
count1 = nob*num_cell;   count2 = count1*Unit_Cells;   N = Unit_Cells;
assert(size(bigH,1) == count1, 'bigH is not 44*num_cell: check ZBNR / proj');
H = hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);

blk = @(s) (s-1)*count1 + (1:count1);
D0 = cell(1,N);  U0 = cell(1,N-1);  L0 = cell(1,N-1);
for s = 1:N
    D0{s} = full(H(blk(s),blk(s)));
    if s < N
        U0{s} = full(H(blk(s),blk(s+1)));       % H_{s,s+1}
        L0{s} = full(H(blk(s+1),blk(s)));       % H_{s+1,s}
    end
end
assert(N < 3 || norm(full(H(blk(1),blk(3))),'fro') < 1e-10, ...
       'H couples non-adjacent slices: RGF needs nearest-slice coupling');

% gate: +lambda on the bottom Te plane (atoms 1,2), -lambda on the top Te plane (atoms 3,4), 0 on W
chi_cell  = [ones(12,1); -ones(12,1); zeros(20,1)];
chi_slice = repmat(chi_cell, num_cell, 1);

if validate_rgf
    z   = (muS+muD)/2 + 1i*eta;   lam = lambda_list(ceil(end/2));
    [g1,~] = lead_eps(bigH,bigUbar,eta,real(z),count1,count2);
    [g2,~] = lead_eps(bigH,bigU,  eta,real(z),count1,count2);
    S1 = bigUbar*g1*bigU;   S2 = bigU*g2*bigUbar;
    T_rgf = rgf_T_trG(z,lam,D0,U0,L0,repmat({chi_slice},1,N),S1,S2,false);
    b1 = blk(1);  bN = blk(N);
    A = z*speye(count2) - sparse(H) - spdiags(lam*repmat(chi_slice,N,1),0,count2,count2);
    A(b1,b1) = A(b1,b1) - S1;   A(bN,bN) = A(bN,bN) - S2;
    RHS = zeros(count2,count1);  RHS(bN,:) = eye(count1);
    X   = A\RHS;                             % last-slice columns of G^r
    G1N = X(b1,:);
    T_full = real(trace(1i*(S1-S1')*G1N*(1i*(S2-S2'))*G1N'));
    fprintf('RGF check at E = %.4f eV, lambda = %.3f eV: T_RGF = %.6f | T_full = %.6f\n', ...
            real(z), lam, T_rgf, T_full);
    T_0 = rgf_T_trG(z,0,D0,U0,L0,repmat({chi_slice},1,N),S1,S2,false);
    fprintf('Pristine, lambda = 0 (device = lead): T = %.6f (must be an integer; 2 in the QSH window)\n', T_0);
    clear A RHS X G1N
end
clear H

%% ===================== Edge vacancies =====================
cell_of  = ceil((1:n_at)'/6);                   % cell of each atom
type_of  = (1:n_at)' - 6*(cell_of-1);           % 1..6 = Te1 Te2 Te3 Te4 W1 W2
slice_of = Device_unit(cell_of,1) + 1;          % slice along x
y_at     = Device_atomic_position(:,2);
y_rows   = unique(y_at);                         % atomic rows across the width
edge_y   = [y_rows(1:n_edge_rows); y_rows(end-n_edge_rows+1:end)];
cand     = find(ismember(y_at,edge_y) & slice_of > n_buffer & slice_of <= N-n_buffer);

rng(seed,'twister');
cand = cand(randperm(numel(cand)));             % one random order -> nested configurations
nV = numel(vac_list);  vac_atoms = cell(1,nV);
atom_name = {'Te1','Te2','Te3','Te4','W1','W2'};
for iv = 1:nV
    Nv = round(vac_list(iv)/100*n_at);
    assert(Nv <= numel(cand), 'Not enough edge sites: increase n_edge_rows');
    vac_atoms{iv} = sort(cand(1:Nv));
    fprintf('vacancy %4s %%: %2d atoms (%.3f %%)', vac_tag{iv}, Nv, 100*Nv/n_at);
    for t = 1:6
        nt = nnz(type_of(vac_atoms{iv}) == t);
        if nt > 0, fprintf('  %s x%d', atom_name{t}, nt); end
    end
    fprintf('\n');
end

% deleted atom: every spin-orbital decoupled and pushed to U_vac (outside the energy window)
Dv = cell(nV,N);  Uv = cell(nV,N-1);  Lv = cell(nV,N-1);  chiv = cell(nV,N);
dg = 1:count1+1:count1^2;
for iv = 1:nV
    rm = cell(1,N);
    for a = vac_atoms{iv}(:)'
        s = slice_of(a);  t = type_of(a);
        rm{s} = [rm{s}, (cell_of(a)-1)*nob - (s-1)*count1 + off(t) + (1:norb(t))];
    end
    for s = 1:N
        r = rm{s};
        D = D0{s};  D(r,:) = 0;  D(:,r) = 0;  D(dg(r)) = U_vac;  Dv{iv,s} = D;
        chi = chi_slice;  chi(r) = 0;  chiv{iv,s} = chi;
        if s < N
            M = U0{s};  M(r,:) = 0;  M(:,rm{s+1}) = 0;  Uv{iv,s} = M;
            M = L0{s};  M(rm{s+1},:) = 0;  M(:,r) = 0;  Lv{iv,s} = M;
        end
    end
end
clear D0 U0 L0

%% ===================== Energy grid =====================
kT_max = kB*max(T_list);
k_bias = round((muS-muD)/dE);
assert(abs(k_bias*dE-(muS-muD)) < 1e-9, 'Vds/dE must be an integer');
k_tail = ceil(n_kT*kT_max/dE);
E  = muD + (-k_tail:k_bias+k_tail)*dE;          % mu_D and mu_S lie on the grid
nE = numel(E);   nL = numel(lambda_list);
assert(E(1) <= Ef-n_kT*kT_max && E(end) >= Ef+n_kT*kT_max, 'grid does not cover the C_Q kernel');
fprintf('%d energies x %d gate points x %d vacancy levels\n', nE, nL, nV);

%% ===================== NEGF sweep: T(E,lambda,vac) and DOS =====================
run_sig = [seed n_edge_rows n_buffer ZBNR Unit_Cells U_vac eta muS Vds];
Trans = zeros(nE,nL,nV);   DOS = zeros(nE,nL,nV);   iE0 = 1;
if resume && exist(ckpt_file,'file') == 2
    S = load(ckpt_file);
    if isequal(S.E,E) && isequal(S.lambda_list,lambda_list) && ...
       isequal(S.vac_list,vac_list) && isequal(S.run_sig,run_sig)
        Trans = S.Trans;  DOS = S.DOS;  iE0 = S.iE_done + 1;
        fprintf('Resuming at energy point %d/%d\n', iE0, nE);
    end
end

t0 = tic;
for iE = iE0:nE                                  % energies are independent
    [g1,~] = lead_eps(bigH,bigUbar,eta,E(iE),count1,count2);    % left lead
    [g2,~] = lead_eps(bigH,bigU,  eta,E(iE),count1,count2);    % right lead
    Sig1 = bigUbar*g1*bigU;
    Sig2 = bigU*g2*bigUbar;
    Tloc = zeros(nL,nV);   Dloc = zeros(nL,nV);
    for iv = 1:nV
        for il = 1:nL
            [Tloc(il,iv), trG] = rgf_T_trG(E(iE)+1i*eta, lambda_list(il), Dv(iv,:), ...
                                   Uv(iv,:), Lv(iv,:), chiv(iv,:), Sig1, Sig2, compute_DOS);
            Dloc(il,iv) = -imag(trG)/(pi*L_dev);     % states/(eV*Angstrom), both spins
        end
    end
    Trans(iE,:,:) = reshape(Tloc,[1 nL nV]);
    DOS(iE,:,:)   = reshape(Dloc,[1 nL nV]);
    if iE == iE0 || mod(iE,10) == 0 || iE == nE
        iE_done = iE;
        save(ckpt_file,'E','lambda_list','vac_list','run_sig','Trans','DOS','iE_done');
        el = toc(t0)/60;
        fprintf('E %3d/%3d (%.4f eV) | %.1f min elapsed | ~%.1f min left\n', ...
                iE, nE, E(iE), el, el*(nE-iE)/(iE-iE0+1));
    end
end

%% ===================== Current and quantum capacitance for each T =====================
Vgs     = lambda_list;                  % same meaning as before (old scripts multiply by 1.26706)
V_TG_BG = gate_factor*lambda_list;      % V
Ez      = lambda_list/t_wte2;           % V/Angstrom (thesis convention)
for iT = 1:numel(T_list)
    Temp = T_list(iT);   kT = kB*Temp;
    for iv = 1:nV
        I = zeros(1,nL);   Do = zeros(1,nL);
        for il = 1:nL
            I(il)  = IE*landauer_window(E, Trans(:,il,iv), muS, muD, kT, n_kT);   % A
            Do(il) = thermal_dos(E, DOS(:,il,iv), Ef, kT, n_kT);             % states/(eV*A)
        end
        Cq2 = q*1e10*Do;                % F/m   (per unit channel length)
        Cq1 = Cq2/(W_rib*1e-10);        % F/m^2 (per unit area)
        vac_percent = vac_list(iv);   vac_atoms_used = vac_atoms{iv};
        tag = sprintf('edge_Vacancy_%s_temp%d', vac_tag{iv}, Temp);
        save(['gate_voltage_current_' tag '.mat'], 'Vgs','V_TG_BG','Ez','I','Length_Device', ...
             'Temp','vac_percent','vac_atoms_used','seed');
        save(['quantum_capacitance_data_' tag '.mat'], 'Vgs','V_TG_BG','Cq1','Cq2','Temp','vac_percent');
        fprintf('T = %3d K | vac = %4s %% | I(0 V) = %.3e A | I(%.2f V) = %.3e A | ratio = %.2e\n', ...
                Temp, vac_tag{iv}, abs(I(1)), V_TG_BG(end), abs(I(end)), abs(I(1))/abs(I(end)));
    end
end
save('edge_vacancy_NEGF_all.mat','E','Trans','DOS','lambda_list','V_TG_BG','Ez','T_list', ...
     'vac_list','vac_tag','vac_atoms','seed','muS','muD','Ef','Vds','L_dev','Length_Device', ...
     'W_rib','n_at','gate_factor','t_ox','eps_ox','n_edge_rows','n_buffer');

%% ===================== Local functions =====================
function [T, trG] = rgf_T_trG(z, lam, D, U, L, chi, Sig1, Sig2, doDOS)
% Recursive Green's function for a block-tridiagonal device.
% D{s} = H_ss, U{s} = H_{s,s+1}, L{s} = H_{s+1,s}; gate on-site = lam*chi{s}.
% T = Tr[Gam1 G_1N Gam2 G_1N'],  trG = Tr[G^r] (only if doDOS).
N = numel(D);   n = size(D{1},1);   Id = eye(n);   dg = 1:n+1:n*n;
gL = cell(1,N);
for s = 1:N
    A = -D{s};
    A(dg) = A(dg) + (z - lam*chi{s}(:)).';          % (E + i*eta) - H_ss - lambda*chi
    if s == 1, A = A - Sig1; else, A = A - L{s-1}*gL{s-1}*U{s-1}; end
    if s == N, A = A - Sig2; end
    gL{s} = A\Id;                                    % left-connected Green's function
    if s == 1, G1N = gL{1}; else, G1N = G1N*U{s-1}*gL{s}; end
end
Gam1 = 1i*(Sig1 - Sig1');   Gam2 = 1i*(Sig2 - Sig2');
T = real(trace(Gam1*G1N*Gam2*G1N'));
trG = 0;
if doDOS
    G = gL{N};   trG = trace(G);
    for s = N-1:-1:1
        G   = gL{s} + gL{s}*U{s}*G*L{s}*gL{s};       % full diagonal block G_ss
        trG = trG + trace(G);
    end
end
end

function J = landauer_window(E, Tr, muS, muD, kT, n_kT)
% J = integral of T(E)[f_S - f_D] dE (eV); T(E) interpolated on a grid finer than kB*T
lo = max(E(1), muD - n_kT*kT);   hi = min(E(end), muS + n_kT*kT);
de = min(kT/10, (E(2)-E(1))/4);
x  = linspace(lo, hi, max(201, ceil((hi-lo)/de)+1));
fS = 1./(1 + exp((x-muS)/kT));   fD = 1./(1 + exp((x-muD)/kT));
J  = trapz(x, interp1(E(:), Tr(:), x(:)).'.*(fS - fD));
end

function Dth = thermal_dos(E, DOS, Ef, kT, n_kT)
% Dth = integral of D(E) F_T(E-Ef) dE, F_T = sech^2((E-Ef)/2kT)/(4kT)  -> C_Q = q^2*Dth
lo = max(E(1), Ef - n_kT*kT);   hi = min(E(end), Ef + n_kT*kT);
de = min(kT/10, (E(2)-E(1))/4);
x  = linspace(lo, hi, max(201, ceil((hi-lo)/de)+1));
FT = sech((x-Ef)/(2*kT)).^2/(4*kT);
Dth = trapz(x, interp1(E(:), DOS(:), x(:)).'.*FT);
end