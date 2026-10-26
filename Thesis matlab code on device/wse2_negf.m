clc
clear 
close all
caption = 'WSe2';

%% //////// Creating the input \\\\\\\\\\\\\\\\\\\\\\
ZNR = 20;  % Zigzag WSe2 nano ribbon, possible values = 4,8,12,16,20,24,28,........
% ZNR counts the number of zigzag chains across the ribbon width
%%%%% increase of znr increases width. 
% Calculating number of cells required
nn = ((ZNR-4)/4)+1;
num_cell = ZNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........
% num_cell refers to how many zigzag rows you add across the width (Y axis)

Ef=0; % Fermi energy calculated using MATLAB code from the WSe2_band.dat file found after wannierization

%% //////// Main Unit cell positions \\\\\\\\\\\\\\\\\\\\\\
xa(1)=0; xa(2)=0; xa(3)=1;       % x-coordinates of the first 3 atoms
ya(1)=0; ya(2)=1; ya(3)=1;       % y-coordinates of the first 3 atoms
for j=0:num_cell-4               % build remaining atoms by shifting pattern
    xa(4+j)=xa(1+j)+1;
    ya(4+j)=ya(1+j)+2;
end
za=zeros(1,length(xa))';         % z = 0 (2D structure)
cp = [xa',ya',za];               % [x,y,z] coordinates for the main unit cell

%% //////// Neighbour Unit cell positions \\\\\\\\\\\\\\\\\\\\\\
xb(1)=1; xb(2)=1; xb(3)=2;       % shifted version of main unit (neighbor)
yb(1)=0; yb(2)=1; yb(3)=1;
for j=0:num_cell-4
    xb(4+j)=xb(1+j)+1;
    yb(4+j)=yb(1+j)+2;
end
zb=zeros(1,length(xb))';
cpp = [xb',yb',zb];              % [x,y,z] for neighbor unit cell%% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\

% \\\\\\\\\\\\\\\\\ Atom position //////////////////
na = num_cell/3;
at1=[1 2 3;0 2 3;1 0 0]; 
at=[1 2 3;0 2 3;1 0 0];
for j=1:na-1
    at=[at;at1];
end

%% \\\\\\\\\\\\\\\ Calling the function///////////////
cell_pos = cp;      % unit cell
cell_pos1 = cpp;    % neighbour unit cell
atomst = at;        % unit cell atoms
atomst1 = at;       % neighbour cell atoms
proj = [10 6 6];    % [n_orb_SeTop, n_orb_W, n_orb_SeBot]
HNN = 20;            % hopping cutoff (nearest-neighbor range)
Vext = 0;           % external voltage = 0 (will add later)

% bigH  = intra-cell Hamiltonian
% bigU  = coupling to next cell
% bigUbar = coupling to previous cell
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);
%[hr,hc]=size(bigH); 



%% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%%% Device calculation
Unit_Cells = 10; %%%% device length 
eta=1e-8;
count1=44*num_cell/3; %%%% number of orbits in a unit cell. 
count2=count1*Unit_Cells; 
Vgs =linspace(-1.2,1,10); % Top or Bottom Gate Voltage;
% Vgs=0.5;
Vds = 0.05;% drain to source voltage Vds
tox = 7; %in angstorm(oxide thickness)
eps_ox = 25; % relative permittivity of HfO2

ts=3.3;eps_s=1.6;ti=3.17;eps_i=3.29;width=4.59; kb=8.3333e-05;Temp=300;% % Temperature
kT = 0.0259;
h_plank = 6.626e-34;
q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05;

IE=(2*q*q/ h_plank); %%%% important 




[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);

%% /////// device calculation \\\\\\\\\\\
Ef = 0; % Fermi energy of Bulk Bi(111) plane 
kT = 0.0259;% at 300 K
eta = 1e-6;
tol = 1e-6;
% I0=(2*q^2)/h;

h = 6.626e-34;

q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h);

I0=(q^2)/hbar/2/pi;

dTe1  =       6.35624338;
dTe2   =      5.73543738;
dTe3    =     9.92929379;
dTe4     =    9.30776171;
dW1      =   7.93340128;
dW2        =  7.73232446;
% % % Temperature
Temp=300;
count=0
% % % Define Fermi function
f=@(E,Temp)1./(1+exp(E/(kb*Temp)));
%% ///////////////// Starting Calculation \\\\\\\\\\\\\\\\\


for kk=1:length(Vgs)

%%%%%% External voltage apply 
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
Voltage=Vgs(kk);

Electric_field=-(2*Voltage)/(dTe2-dTe3);
   Unit_cell_H=zeros(44);


   Ez=Voltage/(2*ti*eps_s/eps_i+ts);
alpha_V=0.675;
lamda_v=alpha_V*Ez;

for i=1:44
    if(i<=6)
    Unit_cell_H(i,i)=-lamda_v;
    end
    if(i<=12&&i>6)
    Unit_cell_H(i,i)=-lamda_v;
    end
    if(i<=18&&i>12)
    Unit_cell_H(i,i)=+lamda_v;
    end
    if(i<=24&&i>18)
    Unit_cell_H(i,i)=+lamda_v;
    end
    
end





BBBB=diag(ones(1,num_cell*Unit_Cells));
big_H_exeternal_voltage=kron(BBBB,Unit_cell_H);
UL=big_H_exeternal_voltage;
U=UL;
Phi=UL;


    mu1=Ef;
    mu2=Ef-Vds;
    E = linspace(mu2,mu1,10);
    dE=E(3)-E(2);
    delphi = 1;
    Ite = 0;
    f1=1./(1+exp((E-mu1)./kT));
	f2=1./(1+exp((E-mu2)./kT));
    epsilon=1;
    rho = 0;

for k1=1:length(E)
    
    [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);

h=bigH;
ud=bigUbar;
u=bigU;
H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
[G1,countG1]=lead_eps(h,ud,eta,E(k1),count1,count2);
[G2,countG2]=lead_eps(h,u,eta,E(k1),count1,count2);
sig1=ud*G1*u;
sig2=u*G2*ud;

% % Calculating and plotting Gn

sig1f=zeros(count2);
sig2f=sig1f;
for m=1:count1
    for n=1:count1
        sig1f(m,n)=sig1(m,n);
        sig2f(m+count2-count1,n+count2-count1)=sig2(m,n);
    end
end


ppp=2
gam1=1i*(sig1f-sig1f');
gam2=1i*(sig2f-sig2f');
sigmain=gam1*f1(k1)+gam2*f2(k1);
sigmaout=gam1*(1-f1(k1))+gam2*(1-f2(k1));
H_total=H+potential_across_length+big_H_exeternal_voltage;
Gr=inv((E(k1)*eye(count2))-H_total-sig1f-sig2f);
Ga=Gr';
Gn=Gr*sigmain*Ga;
    Gp=Ga*sigmaout*Gr;
    rho=rho+(1/(2*pi))*Gn*dE;

   T(k1)=real(trace(gam1*Gr*gam2*Ga));

end
  prho = real(diag(rho));

Electron_density(kk) = sum(diag(rho));
dE=E(2)-E(1);
dI=IE*(f(E-mu1,Temp)-f(E-mu2,Temp)).*T*dE;
I(kk)=sum(dI);
    count=count+1

end


save Igsvsvgs I
hold on
h=plot(Vgs,I*10^6,'b');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[16])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Current Ids (µA) --> ')

