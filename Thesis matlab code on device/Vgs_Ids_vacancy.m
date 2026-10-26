clc
clear 
close all
caption = 'wte2';
%%%%%%%%%% May 2026 , MSc thesis final year ... I aim to finish my MSc.
%%%%%%%%%% Very important code
%%%%%%%%%% May 2026 , I aim to submit a journal paper on device 

%%%%%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 12;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% % % Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% % % //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
Unit_Cells=30;
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
cp = [xa',ya',za];

% % % % % /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
xb(1)=1;xb(2)=1;xb(3)=1;
yb(1)=0;yb(2)=1;yb(3)=2;
for j=0:num_cell-4

    xb(4+j)=xb(1+j)+0;
    yb(4+j)=yb(1+j)+3;
end
zb=zeros(1,length(xb))';
neighbor_unit = [xb', yb', zb];
cpp = [xb',yb',zb];
% 
Device_unit=[];
for i=1:Unit_Cells
   for j = 1:size(main_unit, 1)
    Device_unit(end+1,:)=main_unit(j,:)+[(i-1) 0 0];
   end
end

% 
% 
% % % % % Scaling factors for atoms (example scaling factors for Te and W atoms)
atom_factors = [0.00000044, 1.31473454, 6.35624338;
              1.74897061, 4.43335644, 5.73543738;
             -0.00000002, 5.44760823, 9.92929379;
              1.74897108, 2.22834683, 9.30776171;
              1.74897050, 6.06754902, 7.93340128;
              0.00000040, 3.81415493, 7.73232446];
% 
unit_cell_length=[3.49794  6.33827 15.4319];
ax =  3.49794;
by =  6.33827;
cz =  15.4319;
% 
% % % % % % Initialize atomic positions
main_atomic_position = [];
neighbor_atomic_position=[];
% 
%%%%%%% Generate positions for Te and W atoms
for i = 1:length(main_unit)

       for j = 1:size(atom_factors, 1)
        main_atomic_position(end+1, :) = main_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
       neighbor_atomic_position(end+1, :) = neighbor_unit(i, :) .*unit_cell_length+ atom_factors(j, :);
       end

end
% 
% 
%% %%%%%device atomic position 
Device_atomic_position=[];
for i = 1:length(Device_unit)

       for j = 1:size(atom_factors, 1)
         Device_atomic_position(end+1, :) = Device_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms

       end

end
% 
% %%%% next calculation 
% 


finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
[number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);

x_uniq=unique(finalc_w_extra_lrp(:,1));
y_uniq=unique(finalc_w_extra_lrp(:,2));
z_uniq=unique(finalc_w_extra_lrp(:,3));
xyzlengthuh=[(max(x_uniq)-min(x_uniq)) (max(y_uniq)-min(y_uniq)) (max(z_uniq)-min(z_uniq))]
Vds=0.05;%%%drain voltage
Vgs=linspace(0,1.5,20);%%%%%gate voltage


Length_Device=max(x_uniq)-min(x_uniq);

% 
% % % % % % \\\\\\\\\\\\\\\\\ Atom position //////////////////
na = num_cell/3;
at1=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
at=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
for j=1:na-1
    at=[at;at1];

end

% %%%%%%  \\\\\\\\\\\\\\\ Calling the function///////////////
cell_pos = cp;
cell_pos1 = cpp;
atomst = at;
atomst1 = at;
proj = [6 6 6 6 10 10];
HNN = 5;
Vext = 0;
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 

% 
% % ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
% %%%%%%%%% Device calculation
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells; 
tox = 7; %in angstorm(oxide thickness)
eps_ox = 25; % relative permittivity of HfO2

H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);



%%%%%%%%%%%%%% device vacancy calculation 



% %%%%%%%%%%%% \\\\\\\\\\\\\Edge Vacancy ////////%%%%%%%%%
% 
% 
% % % % Assuming Device_atomic_position is already defined
vector = Device_atomic_position(:,2);  % Extracting the second column

% Find all indices where the value is the minimum
min_val = min(vector);
min_indices = find(vector == min_val);
min_indices1 = find(vector == 2.22834683);

% Find all indices where the value is the maximum
max_val = max(vector);
max_indices = find(vector == max_val);

% Display the results
% disp(['Minimum value is ' num2str(min_val) ' at indices:']);
% disp(min_indices);
% 
% disp(['Maximum value is ' num2str(max_val) ' at indices:']);
% disp(max_indices);

% Choose some of the indices (e.g., randomly selecting 3 from min and max indices)
num_to_select = 5;  % Number of random indices to select

% Randomly select 'num_to_select' indices from min_indices and max_indices
random_min_indices = min_indices(randperm(length(min_indices), num_to_select));
random_max_indices = max_indices(randperm(length(max_indices), num_to_select));

random_min_indices1 = min_indices1(randperm(length(max_indices), 4));

randomNumbers=[random_min_indices; random_max_indices;random_min_indices1];

for i=1:length(randomNumbers)


% Calculate the quotient and remainder
quotient = floor(randomNumbers(i) / 6);  % Integer division
remainder = mod(randomNumbers(i), 6);    % Remainder

if(remainder==1 || remainder==2 || remainder==3 || remainder==4)
    for ii=1:6
        H(quotient*44+ii,:)=0;
        H(:,quotient*44+ii)=0;
    end

end


if(remainder==5 || remainder==0)
    for ii=1:10
        H(quotient*44+ii,:)=0;
        H(:,quotient*44+ii)=0;
    end

end

end
% 
% 
% %%%%% Edge vacancy %%%%%%%









% /////// device calculation \\\\\\\\\\\
Ef=3.015; % Fermi energy of Bulk Bi(111) plane 
kT = 0.0259;% at 300 K
eta = 1e-6;
tol = 0.001;
% 
% 
h = 6.626e-34;

q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h);
I0=(2*q^2)/h;
I0=(q^2)/hbar/2/pi;

dTe1  =       6.35624338;
dTe2   =      5.73543738;
dTe3    =     9.92929379;
dTe4     =    9.30776171;
dW1      =   7.93340128;
dW2        =  7.73232446;
% % % %%%% Temperature
Temp=300;
count=0
% % Define Fermi function
f=@(E,Temp)1./(1+exp(E/(kb*Temp)));
% 
% %%%%%gate voltage
ts=3.3;eps_s=1.6;ti=3.17;eps_i=3.29;width=4.59; kb=8.3333e-05;Temp=300;% % Temperature
kT = 0.0259;
h_plank = 6.626e-34;
q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h_plank);

thbn=3.17;eps_hbn=3.29;twte2=(dTe3-dTe2);eps_wte2=1;

alpha_v=1;

% 
% 
% 
% %% ///////////////// Starting Calculation \\\\\\\\\\\\\\\\\
figure(1)

for kk=1:length(Vgs)

%%%%%% External voltage apply 
% % % % [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);%%%%% [hr,hc]=size(bigH); 
Voltage=Vgs(kk); %%%% the difference between top and bottom gate voltage ... 
Electric_field=(Voltage)/((2*thbn*eps_hbn/eps_wte2)+twte2); 

   Unit_cell_H=zeros(44);

   lamda_v=Voltage;

for i=1:44
    if(i<=6)
    Unit_cell_H(i,i)=+lamda_v;
    end
    if(i<=12&&i>6)
    Unit_cell_H(i,i)=+lamda_v;
    end
    if(i<=18&&i>12)
    Unit_cell_H(i,i)=-lamda_v;
    end
    if(i<=24&&i>18)
    Unit_cell_H(i,i)=-lamda_v;
    end

end


BBBB=diag(ones(1,num_cell*Unit_Cells));
big_H_gate_voltage=kron(BBBB,Unit_cell_H);

for i=1:length(randomNumbers)


% Calculate the quotient and remainder
quotient = floor(randomNumbers(i) / 6);  % Integer division
remainder = mod(randomNumbers(i), 6);    % Remainder

if(remainder==1 || remainder==2 || remainder==3 || remainder==4)
    for ii=1:6
        big_H_gate_voltage(quotient*44+ii,:)=0;
        big_H_gate_voltage(:,quotient*44+ii)=0;
    end

end


if(remainder==5 || remainder==0)
    for ii=1:10
        big_H_gate_voltage(quotient*44+ii,:)=0;
        big_H_gate_voltage(:,quotient*44+ii)=0;
    end

end

end


    mu1=3.05;
    mu2=3.0;
    E = linspace(mu2,mu1,10);
    dE=E(3)-E(2);
    delphi = 1;
    Ite = 0;
    f1=1./(1+exp((E-mu1)./kT));
	f2=1./(1+exp((E-mu2)./kT));
    epsilon=1;
    rho = 0;
    F_T = (1/(4*kT)).*(sech((E-Ef)./(2*kT))).^2;
for k1=1:length(E)
h=bigH;
u=bigU;
ud=bigUbar;
[G1,countG1]=lead_eps(h,ud,eta,E(k1),count1,count2);
[G2,countG2]=lead_eps(h,u,eta,E(k1),count1,count2);
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

gam1=1i*(sig1f-sig1f');
gam2=1i*(sig2f-sig2f');
sigmain=gam1*f1(k1)+gam2*f2(k1);
sigmaout=gam1*(1-f1(k1))+gam2*(1-f2(k1));

Gr=inv((E(k1)*eye(count2))-H-big_H_gate_voltage-sig1f-sig2f);
Ga=Gr';
Gn=Gr*sigmain*Ga;
Gp=Ga*sigmaout*Gr;
T(k1)=trace(gam1*Gr*gam2*Ga);

DOS1(k1) = -(1/(pi*Length_Device))*imag(trace(Gr));
end



dI=IE*(f(E-mu1,Temp)-f(E-mu2,Temp)).*T*dE;
I(kk)=sum(dI);
count=count+1

DOS = DOS1.*F_T;
dE=E(2)-E(1);
Do(kk) = (dE/2)*(DOS(1)+2*sum(DOS(2:end-1))+DOS(end));


end

q = 1.6e-19;
eV = 1.6e-19;
A = 1e-10;
Cq1 = (q^2/(eV*A^2)).*Do; % unit F/m^2
Cq2 = (q^2/(eV*A)).*Do; % unit F/m


save('gate_voltage_current_edge_Vacancy_0.85_temp300.mat','Vgs','I','Length_Device')
save('quantum_capacitance_data_edge_Vacancy_0.85_temp300.mat','Vgs','Cq1','Cq2')
hold on
h=plot(Vgs,I,'b');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Current Ids (A) --> ')

