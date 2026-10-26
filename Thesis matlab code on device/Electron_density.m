clc
clear 
close all
caption = 'wte2';
%%%%%%%%////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 28;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% % Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% % % //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
Unit_Cells=5;

%%%%% zigzag nanoribbon
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

Device_unit=[];
for i=1:Unit_Cells
   for j = 1:size(main_unit, 1)
    Device_unit(end+1,:)=main_unit(j,:)+[(i-1) 0 0];
   end
end



% % % % % Scaling factors for atoms (example scaling factors for Te and W atoms)
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

% % % % % Initialize atomic positions
main_atomic_position = [];
neighbor_atomic_position=[];
% Generate positions for Te and W atoms
for i = 1:length(main_unit)

       for j = 1:size(atom_factors, 1)
        main_atomic_position(end+1, :) = main_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
       neighbor_atomic_position(end+1, :) = neighbor_unit(i, :) .*unit_cell_length+ atom_factors(j, :);
       end

end


%% device atomic position 
Device_atomic_position=[];
for i = 1:length(Device_unit)

       for j = 1:size(atom_factors, 1)
         Device_atomic_position(end+1, :) = Device_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms

       end

end

%%%% next calculation 

Device_atomic_position=Device_atomic_position;

finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
[number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);



% % mod(8,6)
% % % [xxGrid, yyGrid]=meshgrid(x_uniq ,y_uniq);


% % % % % \\\\\\\\\\\\\\\\\ Atom position //////////////////
na = num_cell/3;
at1=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
at=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
for j=1:na-1
    at=[at;at1];

end

% \\\\\\\\\\\\\\\ Calling the function///////////////
cell_pos = cp;
cell_pos1 = cpp;
atomst = at;
atomst1 = at;
proj = [6 6 6 6 10 10];
HNN = 5;
Vext = 0;
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 


% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%% Device calculation
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells; 
tox = 7; %in angstorm(oxide thickness)
eps_ox = 25; % relative permittivity of HfO2



[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);

% /////// device calculation \\\\\\\\\\\
Ef = 3.13195; % Fermi energy of Bulk Bi(111) plane 
kT = 0.0259;% at 300 K
eta = 1e-6;
tol = 0.001;


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
% % Temperature
Temp=300;
count=0
% % Define Fermi function
f=@(E,Temp)1./(1+exp(E/(kb*Temp)));






%% ///////////////// Starting Calculation \\\\\\\\\\\\\\\\\

%% transport calculation 

E2=3.0191;
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells;   
f1=1;
f2=0;
%[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
%H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
h=bigH;
ud=bigUbar;
u=bigU;
H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
[G1,countG1]=lead_eps(h,ud,eta,E2,count1,count2);
[G2,countG2]=lead_eps(h,u,eta,E2,count1,count2);
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
sigmain=gam1*f1+gam2*f2;
sigmaout=gam1*(1-f1)+gam2*(1-f2);

Gr=inv((E2*eye(count2))-H-sig1f-sig2f);
Ga=Gr';
Gn=Gr*sigmain*Ga;
Gp=Ga*sigmaout*Gr;


 T=trace(gam1*Gr*gam2*Ga)


Electron=real(diag(Gn));


 p=2;
  for i=1:(count2/44)
      gn1(i,1)=Electron(p,1);
      p=p+44;

  end

% % 
% % % % save Igsvsvgs I
% % % % hold on
% % % % h=plot(Vgs,I*10^6,'b');
% % % % set(h,'linewidth',[2.0])
% % % % set(gca,'Fontsize',[24])
% % % % grid on
% % % % xlabel(' Voltage Vgs ( V ) --> ')
% % % % ylabel(' Current Ids (µA) --> ')
% % % 
% % % % % 
% % % 
% % % dTe1  =       6.35624338;
% % % dTe2   =      5.73543738;
% % % dTe3    =     9.92929379;
% % % dTe4     =    9.30776171;
% % % dW1      =   7.93340128;
% % % dW2        =  7.73232446;
% % % % % % Extract and scale coordinates
x = finalc_w_extra_lrp(:, 1); % x-coordinates
y = finalc_w_extra_lrp(:, 2); % y-coordinates
z = finalc_w_extra_lrp(:, 3); % z-coordinates


% % % % % % % Define constants and inputs
z_fixed =  6.35624338; % Fixed z-coordinate
z_new = z_fixed; % Optional: Keep z constant or adjust if needed
% Filter for fixed z-coordinate
mask = (z == z_fixed); % Tolerance for fixed z
x_filtered = x(mask); % Filtered x-coordinates
y_filtered = y(mask); % Filtered y-coordinates



% % Optional: Surface fitting with finer grid
% [fitresult, gof] = createFit(x_filtered, y_filtered,gn1 ); % Fit with finer data
% xlabel('Width (x)');
% ylabel('Length (y)');
% zlabel('Potential (V)');
% title('Electron density');
% colorbar;
% axis tight;

LDOS2=Electron;
jjjj=1;
for iiii=1:44:length(LDOS2)
    LDOS(jjjj)=sum(LDOS2(iiii:iiii+6-1));
    LDOS(jjjj+1)=sum(LDOS2(iiii+6:iiii+12-1));
    LDOS(jjjj+2)=sum(LDOS2(iiii+12:iiii+18-1));
    LDOS(jjjj+3)=sum(LDOS2(iiii+18:iiii+24-1));
    LDOS(jjjj+4)=sum(LDOS2(iiii+24:iiii+34-1));
    LDOS(jjjj+5)=sum(LDOS2(iiii+34:iiii+44-1));
    jjjj=jjjj+6;
end
N_cell = Unit_Cells;
count1 = 6*num_cell;

figure()
LDOS1=zeros(N_cell,count1);
kj=0;
for i=1:N_cell
    for j=1:count1
        LDOS1(i,j)=real(LDOS(j+kj));
    end
    kj=kj+count1;
end 


%% 2d plot
LDOS3(1) = 0;
LDOS3(2) = 0;
for i = 1:count1
    LDOS3(i+2)=sum(LDOS1(:,i))/N_cell;
end
LDOS3(end+1) = 0;
LDOS3(end+1) = 0;
x = -1:128;
[fitresult, gof] = createFit(x, LDOS3);


save zigzag LDOS3

