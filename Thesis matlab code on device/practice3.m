clc 
close all 
clear all
caption = 'wte2';
%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 4;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required

nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........



% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
%% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%%% Device calculation
Unit_Cells = 5; %Number of Unit Cells
eta=1e-8;
pi=3.1416;
count1=44*num_cell;
count2=count1*Unit_Cells; 

unit_cells=Unit_Cells;
ax =  3.49794;
by =  6.33827;
cz =  15.4319;


xa(1)=0;xa(2)=0;xa(3)=0;
ya(1)=0;ya(2)=1;ya(3)=2;
for j=0:num_cell-4
    
    xa(4+j)=xa(1+j)+0;
    ya(4+j)=ya(1+j)+3;
end
za=zeros(1,length(xa))';

main_unit = [xa', ya', za];

% /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
xb(1)=1;xb(2)=1;xb(3)=1;
yb(1)=0;yb(2)=1;yb(3)=2;
for j=0:num_cell-4
    
    xb(4+j)=xb(1+j)+0;
    yb(4+j)=yb(1+j)+3;
end
zb=zeros(1,length(xb))';

neighbor_unit = [xb', yb', zb];


Device_unit=[];
for i=1:unit_cells
   for j = 1:size(main_unit, 1)
    Device_unit(end+1,:)=main_unit(j,:)+[(i-1) 0 0];
   end
end



% Scaling factors for atoms (example scaling factors for Te and W atoms)
atom_factors = [0.00000044, 1.31473454, 6.35624338;
              1.74897061, 4.43335644, 5.73543738;
             -0.00000002, 5.44760823, 9.92929379;
              1.74897108, 2.22834683, 9.30776171;
              1.74897050, 6.06754902, 7.93340128;
              0.00000040, 3.81415493, 7.73232446];

unit_cell_length=[3.49794  6.33827 15.4319];
ax =  3.49794;
by =  6.33827;
cz =  15.4319;

% Initialize atomic positions
main_atomic_position = [];
neighbor_atomic_position=[];
% Generate positions for Te and W atoms
for i = 1:length(main_unit)

       for j = 1:size(atom_factors, 1)
        main_atomic_position(end+1, :) = main_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
       neighbor_atomic_position(end+1, :) = neighbor_unit(i, :) .*unit_cell_length+ atom_factors(j, :);
       end
    
end


%%% device atomic position 
Device_atomic_position=[];
for i = 1:length(Device_unit)

       for j = 1:size(atom_factors, 1)
         Device_atomic_position(end+1, :) = Device_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
     
       end
    
end



%%%%% next calculation 

Device_atomic_position=Device_atomic_position.*10^-10;

finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
[number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);



caption = 'wte2';
%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 4;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
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
%[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 


%% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%%% Device calculation
Unit_Cells = 5; %Number of Unit Cells
eta=1e-8;
pi=3.1416;
count1=44*num_cell;
count2=count1*Unit_Cells; 
% Vgs =linspace(0,1.1,10); % Top or Bottom Gate Voltage;
Vgs=[0.1];
n0=0;
Vds = 0.2;% drain to source voltage Vds
% Material-specific constants for WTe2
a = 3.5e-10; % lattice constant in meters
K = 10; % dielectric constant (placeholder value, needs actual value)
r = 5e-9; % radius in meters
tox = 5e-9; % oxide thickness in meters
q = 1.6e-19; % elementary charge in Coulombs
eps0 = 8.854e-12; % permittivity of free space in F/m

% Calculate initial potential U0 for WTe2
U0 = q / (2 * pi * a * K * eps0) * log((r + tox) / r);




[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);







%% /////// device calculation \\\\\\\\\\\
Ef = 3.4055; % Fermi energy of Bulk Bi(111) plane 
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
% Temperature
Temp=300;
count=0
% Define Fermi function
f=@(E,Temp)1./(1+exp(E/(kb*Temp)));
%% ///////////////// Starting Calculation \\\\\\\\\\\\\\\\\


for kk=1:length(Vgs)

%%%%%% External voltage apply 
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
Voltage=Vgs(kk);

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


BBBB=diag(ones(1,num_cell*Unit_Cells));
big_H_exeternal_voltage=kron(BBBB,Unit_cell_H);
UL=diag(big_H_exeternal_voltage);
U=UL;
Phi=UL;


    mu1=Ef;
    mu2=Ef-Vds;
    E = linspace(mu2-4*kT,mu1+4*kT,20);
    dE=E(3)-E(2);
    delphi = 1;
    Ite = 0;
    f1=1./(1+exp((E-Ef)./kT));
	f2=1./(1+exp((E-Ef+Vds)./kT));
    epsilon=1;
    rho = 0;
    rho = 0;
    Phi_old = Phi;

for k1=1:length(E)
%     f1=1/(1+exp((E(k1)-mu1)/kT));f2=1/(1+exp((E(k1)-mu2)/kT));
    [G1,countG1]=lead_eps(bigH,bigUbar,eta,E(k1),count1,count2);
    [G2,countG2]=lead_eps(bigH,bigU,eta,E(k1),count1,count2);
    sig1=bigUbar*G1*bigU;
    sig2=bigU*G2*bigUbar;
    sig1f=zeros(count2);
    sig2f=sig1f;
    sig1f(1:count1,1:count1)=sig1;
    sig2f(count2-count1+1:end,count2-count1+1:end)=sig2;
    gam1=1i*(sig1f-sig1f');
    gam2=1i*(sig2f-sig2f');
    sigmain=gam1*f1(k1)+gam2*f2(k1);
    sigin1=f1(k1)*gam1;sigin2=f2(k1)*gam2;
    Gr=inv((E(k1)*eye(count2))-H-sig1f-sig2f-diag(Phi));A=i*(Gr-Gr');
    Ga=Gr';
    Gn=Gr*sigmain*Ga;
    rho=rho+(1/(2*pi))*Gn*dE;

   T(k1)=real(trace(gam1*Gr*gam2*Ga));
   I1(k1)=real(trace(sigin1*A)-trace(gam1*Gn));
  I2(k1)=-real(trace(sigin2*A)-trace(gam2*Gn));
  rho=rho+(1/(2*pi))*Gn*dE;

end
  prho = real(diag(rho));
end




%% ////////// Adding extra points \\\\\\\\\\\\\\\\\\\\\\\\
 k11=1;
 x_max=max(Device_atomic_position(:,1));
x_min=min(Device_atomic_position(:,1));
 % % % ch_r = linspace(x_min,x_max,40);

%%%  % ch_r = 5.73543738:dz:9.92929379;
 ch_x=linspace(x_min*10^-10,x_max*10^-10,40);
 dx=ch_x(2)-ch_x(1);
% 
 n_extra = length(ch_x);
for i=1:length(finalc_w_extra_lrp)

    for j=0:n_extra-1
    finalc_x(k11+j,:)=finalc_w_extra_lrp(i,:);
    end
    k11=k11+n_extra;
end

for i=1:n_extra:length(finalc_x)
    for j=1:n_extra
        finalc_x(i+j-1,1)=ch_x(j);
    end

end
% 
% 
% 



%%%% add extra point in y direction 

%% ////////// Adding extra points \\\\\\\\\\\\\\\\\\\\\\\\
 k11=1;
 y_max=max(Device_atomic_position(:,2));
y_min=min(Device_atomic_position(:,2));
 % % % ch_r = linspace(x_min,x_max,40);

%%%  % ch_r = 5.73543738:dz:9.92929379;
 ch_y=linspace(y_min*10^-10,y_max*10^-10,40);
 dy=ch_y(2)-ch_y(1);
% 
 n_extra = length(ch_y);
for i=1:length(finalc_w_extra_lrp)

    for j=0:n_extra-1
    finalc_y(k11+j,:)=finalc_w_extra_lrp(i,:);
    end
    k11=k11+n_extra;
end

for i=1:n_extra:length(finalc_y)
    for j=1:n_extra
        finalc_y(i+j-1,1)=ch_y(j);
    end

end
% 
% % 
% % 
% % 
% % %% Calling the system differential matrix function 
% % 
% % 
% %  every wte2 has [6 6 6 6 10 10] orbital
j=1;
for i=1:44:length(prho)
    prho_new(j)=-sum(prho(i:i+6-1));
    prho_new(j+1)=-sum(prho(i+6:i+12-1));
    prho_new(j+2)=-sum(prho(i+12:i+18-1));
    prho_new(j+3)=-sum(prho(i+18:i+24-1));
    prho_new(j+4)=-sum(prho(i+24:i+34-1));
    prho_new(j+5)=-sum(prho(i+34:i+44-1));
    j=j+6;
end
A_X1=Device_atomic_position(:,1);
A_Y1=Device_atomic_position(:,2);
A_Z1=Device_atomic_position(:,3);
f1D_x = griddata(A_X1(:,1),A_Y1(:,1),A_Z1(:,1),prho_new',finalc_x(:,1),finalc_x(:,2),finalc_x(:,3),'nearest');
f1D_y = griddata(A_X1(:,1),A_Y1(:,1),A_Z1(:,1),prho_new',finalc_y(:,1),finalc_y(:,2),finalc_y(:,3),'nearest');

x=ch_x;
y=ch_y;

[X, Y] = meshgrid(ch_x, ch_y); % 2D meshgrid
Xu = X(:); % Unravel 2D meshgrid to 1D array
Yu = Y(:);

% Dirichlet boundary conditions
uL = Vds; % Left
uR = 0;   % Right
uT = Vgs;  % Top
uB = -Vgs;  % Bottom

% Define independent variables
Nx = length(ch_x);  % Number of grid points in x direction
Ny = length(ch_y);  % Number of grid points in y direction
% Boundary indices
ind_L = find(Xu == x(1));    % Left boundary
ind_R = find(Xu == x(end));  % Right boundary
ind_B = find(Yu == y(1));    % Bottom boundary
ind_T = find(Yu == y(end));  % Top boundary



% Source function (right-hand side vector)
g = zeros(Nx * Ny, 1);

% Finite difference operators
Dx_2d = spdiags([-ones(Nx,1), ones(Nx,1)] / (2*dx), [-1, 1], Nx, Nx);
Dy_2d = spdiags([-ones(Ny,1), ones(Ny,1)] / (2*dy), [-1, 1], Ny, Ny);
D2x_2d = spdiags([ones(Nx,1), -2*ones(Nx,1), ones(Nx,1)] / dx^2, [-1, 0, 1], Nx, Nx);
D2y_2d = spdiags([ones(Ny,1), -2*ones(Ny,1), ones(Ny,1)] / dy^2, [-1, 0, 1], Ny, Ny);

% 2D Laplacian operator
Ix = speye(Nx);
Iy = speye(Ny);
L_sys = kron(Iy, D2x_2d) + kron(D2y_2d, Ix);
% Construct right-hand side vector
b = g;
b(ind_L) = uL;
b(ind_R) = uR;
b(ind_T) = uT;
b(ind_B) = uB;

% Solve the sparse system
u = reshape(L_sys \ b, Ny, Nx);



% eps0 = 8.854e-12;
% eps = 15*eps0;
% f1D = (f1D.*1e-19)./(eps); % Check again
% 
% 
% 
% %% Calling Differential matrix 
% [H_poison] = Diffmat_1D( n_extra, number_atoms,dz);
% 
% % 
% % 
% % 
% for ii=1:length(f1D)
%      if(finalc(ii,3)==x_max*10^-10)
%          f1D(ii)=vds;
%      end
%      if(finalc(ii,3)==x_min*10^-10)
%          f1D(ii)=0;
%      end
% 
% end
% % 
% potential =inv(H_poison)*(f1D);
% % 
% for ss=1:(length(potential)/n_extra)
%     distance_etter=100;
%     queryX = A_X1(ss, 1);
%     queryY = A_Y1(ss, 1);
%     queryZ = A_Z1(ss, 1);
% 
%     start=(ss-1)*n_extra+1;
% 
%     for kk=start:start+n_extra-1
%     distances_cal = sqrt((finalc(kk,1) - queryX)^2 + (finalc(kk,2) - queryY)^2 + (finalc(kk,3) - queryZ)^2);
%     if(distances_cal<distance_etter)
%     potential_new(ss)=potential(kk);
%     end
% 
%     distance_etter=distances_cal;
%     end
% end
% 
% hhhh=0;
% 
% % 
% for ll=1:6:length(potential_new)
% 
%       for gg=1:6
%       potential_new_input(gg+hhhh,1)=potential_new(ll);
%       end
%       for gg=7:12
%       potential_new_input(gg+hhhh,1)=potential_new(ll+1);
%       end
%       for gg=13:18
%       potential_new_input(gg+hhhh,1)=potential_new(ll+2);
%       end
%       for gg=19:24
%       potential_new_input(gg+hhhh,1)=potential_new(ll+3);
%       end
%       for gg=25:34
%       potential_new_input(gg+hhhh,1)=potential_new(ll+4);
%       end
%       for gg=35:44
%       potential_new_input(gg+hhhh,1)=potential_new(ll+5);
%       end
%       hhhh=hhhh+44;
% end
% % 
% % 
% % 
