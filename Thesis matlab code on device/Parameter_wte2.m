clc
clear 
close all
caption = 'wte2';
%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
Zigzag_atoms = 16;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required
cal=0
Unit_cell_underlap=0


ZBNR=Zigzag_atoms;
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
Underlap_source=Unit_cell_underlap;
Underlap_drain=Unit_cell_underlap;
Channel_unit_cell=15;

Unit_Cells=Underlap_source+Channel_unit_cell+Underlap_drain;
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

% /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
xb(1)=1;xb(2)=1;xb(3)=1;
yb(1)=0;yb(2)=1;yb(3)=2;
for j=0:num_cell-4

    xb(4+j)=xb(1+j)+0;
    yb(4+j)=yb(1+j)+3;
end
zb=zeros(1,length(xb))';
neighbor_unit = [xb', yb', zb];
cpp = [xb',yb',zb];

Device_unit=[];
for i=1:Unit_Cells
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
% % Generate positions for Te and W atoms
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

%%%% next calculation 

Device_atomic_position=Device_atomic_position;

finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
[number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);

x_uniq=unique(finalc_w_extra_lrp(:,1));
y_uniq=unique(finalc_w_extra_lrp(:,2));
z_uniq=unique(finalc_w_extra_lrp(:,3));
xyzlengthuh=[(max(x_uniq)-min(x_uniq)) (max(y_uniq)-min(y_uniq)) (max(z_uniq)-min(z_uniq))]
Vds=0.05;%%%drain voltage

Unit_cell_bias_drain=zeros(44);
Unit_cell_bias_source=zeros(44);
for i=1:44
    Unit_cell_bias_drain(i,i)=Vds;
    Unit_cell_bias_source(i,i)=0;

end

BBBB=diag(ones(1,num_cell));
big_H_bias_drain=kron(BBBB,Unit_cell_bias_drain);
big_H_bias_source=kron(BBBB,Unit_cell_bias_source);

Length_Device=max(x_uniq)-min(x_uniq);
Width_Device=max(y_uniq)-min(y_uniq);

potential_across_length=zeros(44*num_cell*Unit_Cells,44*num_cell*Unit_Cells);
potential_across_hight=zeros(number_atoms,1);
oooo=0;
for iii=1:number_atoms
       hi=Vds*(finalc_w_extra_lrp(iii,1)-min(x_uniq))/Length_Device;
       if(mod(iii,6)==1||mod(iii,6)==4)
       for jjj=1:6
           potential_across_length(jjj+oooo,jjj+oooo)=hi;

       end
       oooo=oooo+6;
       end

       if(mod(iii,6)==2||mod(iii,6)==3)
       for jjj=1:6
           potential_across_length(jjj+oooo,jjj+oooo)=hi;

       end
       oooo=oooo+6;
       end

       if(mod(iii,6)==5||mod(iii,6)==0)
       for jjj=1:10
           potential_across_length(jjj+oooo,jjj+oooo)=hi;
       end
       oooo=oooo+10;
       end
end






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


%% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%%% Device calculation
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells; 
Vgs =linspace(-12,0,40); % Top or Bottom Gate Voltage;
% Vgs=0.5;
Vds = 0.05;% drain to source voltage Vds
tox = 7; %in angstorm(oxide thickness)
eps_ox = 25; % relative permittivity of HfO2

ts=4.193856;eps_s=1.6;ti=3.17;eps_i=3.29;width=4.59; kb=8.3333e-05;Temp=300;% % Temperature
kT = 0.0259;
h_plank = 6.626e-34;
q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h_plank);




[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);

%% /////// device calculation \\\\\\\\\\\
Ef = 3.13195; % Fermi energy of Bulk Bi(111) plane 
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

for i=1:Underlap_source*num_cell
    BBBB(i,i)=0;
end


for i=1:Underlap_drain*num_cell
    BBBB(end-i,end-i)=0;
end


big_H_exeternal_voltage=kron(BBBB,Unit_cell_H);
UL=big_H_exeternal_voltage;
U=UL;
Phi=UL;

    E = Ef;
    
    f1=1;
	f2=0;
    epsilon=1;
    rho = 0;

    [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);

h=bigH;
ud=bigUbar;
u=bigU;
H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
[G1,countG1]=lead_eps(h,ud,eta,E,count1,count2);
[G2,countG2]=lead_eps(h,u,eta,E,count1,count2);
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


ppp=2;
gam1=1i*(sig1f-sig1f');
gam2=1i*(sig2f-sig2f');
sigmain=gam1*f1+gam2*f2;
sigmaout=gam1*(1-f1)+gam2*(1-f2);
H_total=H+potential_across_length+big_H_exeternal_voltage;
Gr=inv((E*eye(count2))-H_total-sig1f-sig2f);
Ga=Gr';
Gn=Gr*sigmain*Ga;
    Gp=Ga*sigmaout*Gr;

   T(kk)=real(trace(gam1*Gr*gam2*Ga));
q = 1.6e-19; % %%%%Elementary charge (C)
DOS(kk)=-1/pi*imag(trace(Gr))/q/(xyzlengthuh(1)*xyzlengthuh(2)*1e-10*1e-10);
Cq(kk)=DOS(kk)*q^2*1e15*1e-12;
  

end




% % %%%calculation of transcondutance 
gm=IE*Vds*(diff(T)./diff(Vgs))/Length_Device*(1e3*1e3);%%%%% mS/um

figure(1)
h = plot(Vgs(1:end-1), gm, 'b'); 
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[16])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Transconductance(mS/\mum) --> ')

%%% cutoff frequency
% % % (fF/\mum^2)
%% Cox
% Parameters
epsilon_0 = 8.854e-12; % Vacuum permittivity (F/m)
epsilon_r = 25; % Relative permittivity
tox = 3.17e-9; % Oxide thickness (m)

% Oxide capacitance calculation
Cox = (epsilon_0 * epsilon_r / tox)*1e15*1e-12; % Oxide capacitance per unit area (uF/cm^2)

Cgs= (Cox.^-1+Cq.^-1+Cox.^-1).^-1;
Cgd=Cgs(1:end-1);
cut_freq=abs(gm./(2*pi*(Cgd+Cgd)))*(1e-3*1e15);

FT=mean(cut_freq)
figure(2)
h = plot(Vgs(1:end-1),cut_freq, 'b'); 
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[16])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Cut of frequency (HZ) ')


Delay_time=1./(2*pi*cut_freq)
figure(3)
h = plot(Vgs(1:end-1),Delay_time, 'b'); 
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[16])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Delay time (s) ')

Delay=mean(Delay_time);



%% Plot Total Capacitance
figure(4);  % Open figure window 7

% Plot the various capacitances
plot(Vgs, Cox*ones(length(Vgs),1), 'r', 'LineWidth', 2);  % Cox
hold on;  % Hold for subsequent plots
plot(Vgs, Cq, 'g--', 'LineWidth', 2);  % Cq
plot(Vgs, Cg, 'b--', 'LineWidth', 2);  % Cg

% Add labels and title
xlabel('V_{GS} (V)');
ylabel('Capacitances (fF/\mum^2)');
title('Capacitances vs V_{GS} for NR TIFET');
% Add a legend for the three plots
legend('C_{ox}','C_q','C_g', 'Location', 'best');  % Improved labels

% Turn on grid and set axes limits if needed
grid on;


