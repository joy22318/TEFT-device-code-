clc;clearvars;close all;

%% Parameters
e = 1.6e-19; % Electron charge (C)
Vth = 0.2; % Threshold voltage (V)
Vds = 0.1; % Drain-source voltage (V)
Vgs = 0:0.2:1.1; % Gate-source voltage
dVgs=0.01;
f1=1*e;
f2=0*e;
mu = 0.1; % Chemical potential (eV)
h = 6.67e-34; % Reduced Planck's constant (J*s)
a = 1e-9; % Lattice constant (m)
N = 432; % Number of sites
V_x= zeros(N, N);
I=eye(N);
to=8.20*e;
kT=2.59e-3*e;
beta=1/kT;
Wi=2.6e-9*1e6;
L=7.8e-9*1e6;
hbar = 1.0545718e-34; % Reduced Planck's constant in J*s

%% Hamiltonian for a simple 1D chain
%data= load('hamiltonian.mat');
%H= data.bigH;

%% Calculate Ids
Ids = zeros(size(Vgs));
Cq = zeros(size(Vgs));
C = zeros(size(Vgs));
T = zeros(size(Vgs));
m_eff = zeros(size(Vgs));

caption = 'wte2';

%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 12;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
xa(1)=0;xa(2)=0;xa(3)=0;
ya(1)=0;ya(2)=1;ya(3)=2;
for j=0:num_cell-4

    xa(4+j)=xa(1+j)+0;
    ya(4+j)=ya(1+j)+3;
end
za=zeros(1,length(xa))';
cp = [xa',ya',za];

% /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
xb(1)=1;xb(2)=1;xb(3)=1;
yb(1)=0;yb(2)=1;yb(3)=2;
for j=0:num_cell-4

    xb(4+j)=xb(1+j)+0;
    yb(4+j)=yb(1+j)+3;
end
zb=zeros(1,length(xb))';
cpp = [xb',yb',zb];

% \\\\\\\\\\\\\\\\\ Atom position //////////////////
na = num_cell/3;
at1=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
at=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
for j=1:na-1
    at=[at;at1];

end

%% \\\\\\\\\\\\\\\ Calling the function///////////////
cell_pos = cp;
cell_pos1 = cpp;
atomst = at;
atomst1 = at;
proj = [6 6 6 6 10 10];
HNN = 5;
Vext = 0;
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 

% %%
% kpath=[0 0.0 0;
%      0 0.5 0;
%      0.0 0.0 0];
kpath=[0 -0.5 0;
     0 0.0 0;
     0 0.5 0];

% kpoint_fn=[1 20 39];
nopot=[30 30];

% kpath=[0.000000000    0.000000000    0.000000000;
%        0.500000000    0.000000000    0.000000000;
%        0.666666667   -0.333333333    0.000000000;
%        0.000000000    0.000000000    0.000000000];
% kpoint_fn=[1 20 39];
% nopot=[20 10 20];

[kpoints,kpo]=kpgen(kpath,nopot);
count=0
 % Top or Bottom Gate Voltage;

q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(q*q)/(2*pi*hbar);
I0=(q^2)/hbar/2/pi;
dTe1  =       6.35624338;
dTe2   =      5.73543738;
dTe3    =     9.92929379;
dTe4     =    9.30776171;
dW1      =   7.93340128;
dW2        =  7.73232446;


for ser = 1:length(Vgs)

%%%%%% External voltage apply 
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
Voltage=Vgs(ser);


Electric_field=-(2*Voltage)/(dTe2-dTe3);
   Unit_cell_H=zeros(44);
for i=1:44
    if(i<=6)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dTe1)+Voltage;
    end
    if(i<=12&&i>6)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dTe2)+Voltage;
    end
    if(i<=18&&i>12)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dTe3)+Voltage;
    end
    if(i<=24&&i>18)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dTe4)+Voltage;
    end
    if(i<=34&&i>24)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dW1)+Voltage;
    end
    if(i<=44&&i>34)
    Unit_cell_H(i,i)=Electric_field*(dTe2-dW2)+Voltage;
    end
end


BBBB=diag(ones(1,num_cell));
big_H_exeternal_voltage=kron(BBBB,Unit_cell_H);

bigH=bigH+big_H_exeternal_voltage;

%%%% External voltage apply end 

        HB=cell(3,1);
        HB{1,1}=sparse(bigUbar);
        HB{2,1}=sparse(bigH);
        HB{3,1}=sparse(bigU);
        unique_trans_vec=[0 -1 0;0 0 0;0 1 0];
        [hr,hc]=size(HB{2,1});
        diag(HB{2,1});

        Es=[];
        [kpr,kpc]=size(kpoints);
        xx=1;
        yy=1;
        for k=1:kpr
            Hk=zeros(hr,hc);
            for j=1:3
                Hb=HB{j,1};
                Hk=Hk+HB{j,1}.*exp(1i*sum(unique_trans_vec(j,:).*kpoints(k,:)));
            end
            [W,Ek]=eig(Hk);
                 En=sort(real(diag(Ek)));
                Es(:,k)=En;                       
        end

        Ess=Es';


q = 1.6e-19; % Charge of an electron in Coulombs     
b = 169; % Band index of interest
dim = 2; % k-direction of interest (for example, k_x, k_y, or k_z)

% First derivative dE/dk
dE_dk = diff(Ess(:, b)) ./ diff(kpoints(:, dim));

% Second derivative d^2E/dk^2
dE_dk2 = diff(dE_dk) ./ diff(kpoints(2:end, dim));
    m_eff(ser) = hbar^2 / mean(dE_dk2)/q^2; % Effective mass in kg 
    k_B = 1.38064852*(10^(-23)); % the boltzmann constant, Unit: m2*kg*s-2 K-1
    Temp = 300; % the temperature, Unit: K
    DOS2 = (m_eff(ser)/(2*pi*hbar^2));
    C(ser)=DOS2*q^2*1e15*1e-12;

    % Calculate Density of States (DOS) for nanowire
    E_range = linspace(min(Ess(:)), max(Ess(:)), 1000); % Energy range for DOS
    m_eff_a=2.3507e-30;
    DOS = (m_eff_a/(2*pi*hbar^2));
    E_G = 0.366*q; % J = kg*m2*s-2
    Cq = e^2*DOS*(1+exp(E_G/(2*k_B*Temp)./(2*cosh(e*Vgs./(k_B*Temp))))).^(-1);


%% transport calculation 
Unit_Cells=4;
E2=3.39;
Ef=3.4055;
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells;  
mu1=Ef;
mu2=Ef-Vds;
f1=0;
f2=Vds;
%[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
%H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
hd=bigH;
ud=bigUbar;
u=bigU;
H= hamiltonian_TMD_device_H(hd,u,Unit_Cells);
[G1,countG1]=lead_eps(hd,ud,eta,E2,count1,count2);
[G2,countG2]=lead_eps(hd,u,eta,E2,count1,count2);
sig1=ud*G1*u;
sig2=u*G2*ud;

sig1f=zeros(count2);
sig2f=sig1f;
for m=1:count1
    for n=1:count1
        sig1f(m,n)=sig1(m,n);
        sig2f(m+count2-count1,n+count2-count1)=sig2(m,n);
    end
end

ppp=2;
gam1=1i*(sig1f-sig1f');
gam2=1i*(sig2f-sig2f');
sigmain=gam1*f1+gam2*f2;
sigmaout=gam1*(1-f1)+gam2*(1-f2);

 Gr=inv((E2*eye(count2))-H-sig1f-sig2f);
 Ga=Gr';
 Gn=Gr*sigmain*Ga;
 Gp=Ga*sigmaout*Gr;
 transmittence=trace(gam1*Gr*gam2*Ga);
 T(ser)=0.5*abs(transmittence);
 Ids(ser) = (2 * e*e / h) * T(ser) * Vds*1e6; % Calculate current using Landauer formula
 count=count+1
end

%% Determine I_OFF & I_ON
% I_ON = Ids(21);
% I_OFF = Ids(45);

%% Plot Ids vs Vgs
figure(1);
%plot(Vgs+0.4, (Vgs>-0.4).*abs(Ids)+(Vgs<=-0.4).*0.957157, 'LineWidth', 2);
plot(Vgs,Ids,'LineWidth', 2);
hold on;
% plot(0.24, I_OFF, 'ro', 'MarkerSize', 8, 'DisplayName', 'I_{OFF}');
% plot(0, I_ON, 'go', 'MarkerSize', 8, 'DisplayName', 'I_{ON}');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (uA)');
title('I_{DS} vs. V_{GS} for Topological Insulator Channel Nano MOSFET');
% legend('I_{DS}', 'I_{OFF}', 'I_{ON}');
grid on;
%xlim([-0.09 0.24]);
%ylim([0 1]);
%figure(2);plot(Vgs,T,'LineWidth', 2);


%% Plot Cq vs Vgs
figure(2);
plot(Vgs, C, 'LineWidth', 2);
xlabel('V_{GS} (V)');
ylabel('Cq (fF/um^2)');
title('Cq vs. V_{GS} for Topological Insulator Channel Nano MOSFET');
legend('Cq');
grid on;

%% Cox
% Parameters
epsilon_0 = 8.854e-12; % Vacuum permittivity (F/m)
epsilon_r = 2.2; % Relative permittivity
tox = 2e-9; % Oxide thickness (m)

% Oxide capacitance calculation
Cox = (epsilon_0 * epsilon_r / tox)*1e15*1e-12; % Oxide capacitance per unit area (uF/cm^2)

% Plotting
figure(3);
plot(Vgs, Cox * ones(size(Vgs)), 'LineWidth', 2);
xlabel('V_{GS} (V)');
ylabel('C_{ox} (fF/um^2)');
title('Oxide Capacitance vs. V_{GS} for SiO2 Oxide Layer');
grid on;

%% Cg

Cg= (C.^-1+Cox.^-1++Cox.^-1).^-1;
figure(4);
plot(Vgs, Cg, 'LineWidth', 2);
xlabel('V_{GS} (V)');
ylabel('C_{g} (fF/um^2)');
title('Gate Capacitance vs. V_{GS} for Porous SiO2 Oxide Layer');
grid on;

%% Cut_off frequency

% Tr=[T_bar T1(1:41)];
%gm= (2*e*e/h).*Tr*1e-3*(1/dVgs);
%gm= (2*e*e/h)*0.5.*T*Vds;
rho=0.05;    %0.05 ohm-cm at 7*10^18 doping [usable]
Lg= 10*3.4979*1e-7;    % 5nm
Wg= 6*6.33827*1e-7;    % 3nm
tg= 2e-7;    % 2nm
Rg= (rho*Lg)/(Wg*tg);
tau= Rg.*Cg*Wi*L;
f_cut= (1./(2*pi*tau))*1e-9*1e15;

figure(5);
hold on;

plot(Vgs, f_cut, 'LineWidth', 2);

xlabel('V_{GS} (V)');
ylabel('Cutoff Frequency, f_{cut} (GHz)');
title('f_{cut} vs V_{gs} for NR TIFET ');
grid on;

delay=tau*1e-15*1e12;

figure(6);
hold on;

plot(Vgs, delay, 'LineWidth', 2);

xlabel('V_{GS} (V)');
ylabel('Delay (pS)');
title('Delay vs V_{gs} for NR TIFET ');
grid on;

%% Plot Total Capacitance
figure(7);  % Open figure window 7

% Plot the various capacitances
plot(Vgs, Cox*ones(length(Vgs),1), 'r', 'LineWidth', 2);  % Cox
hold on;  % Hold for subsequent plots
plot(Vgs, C, 'g--', 'LineWidth', 2);  % Cq
plot(Vgs, Cg, 'b--', 'LineWidth', 2);  % Cg

% Add labels and title
xlabel('V_{GS} (V)');
ylabel('Capacitances (fF/\mum^2)');
title('Capacitances vs V_{GS} for NR TIFET');

% Add a legend for the three plots
legend('C_{ox}','C_q','C_g', 'Location', 'best');  % Improved labels

% Turn on grid and set axes limits if needed
grid on;

% Optionally, you can uncomment the line below to set x-axis limits
% xlim([-0.6 0.2]);

% Calculate the average effective mass if needed
m_eff_avg = sum(m_eff) / 59;



mean(delay)
mean(f_cut)
mean(Cg)
