% clc
% clear 
% close all
% caption = 'wte2';
% % ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
% ZBNR = 8;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% % % % Calculating number of cells required
% nn = ((ZBNR-4)/4)+1;
% num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........
% 
% % % % //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
% Unit_Cells=6;
% ax =  3.49794;
% by =  6.33827;
% cz =  15.4319;
% xa(1)=0;xa(2)=0;xa(3)=0;
% ya(1)=0;ya(2)=1;ya(3)=2;
% for j=0:num_cell-4
% 
%     xa(4+j)=xa(1+j)+0;
%     ya(4+j)=ya(1+j)+3;
% end
% za=zeros(1,length(xa))';
% main_unit = [xa', ya', za];
% cp = [xa',ya',za];
% 
% % % % % % /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
% xb(1)=1;xb(2)=1;xb(3)=1;
% yb(1)=0;yb(2)=1;yb(3)=2;
% for j=0:num_cell-4
% 
%     xb(4+j)=xb(1+j)+0;
%     yb(4+j)=yb(1+j)+3;
% end
% zb=zeros(1,length(xb))';
% neighbor_unit = [xb', yb', zb];
% cpp = [xb',yb',zb];
% 
% Device_unit=[];
% for i=1:Unit_Cells
%    for j = 1:size(main_unit, 1)
%     Device_unit(end+1,:)=main_unit(j,:)+[(i-1) 0 0];
%    end
% end
% 
% 
% 
% % % % % % Scaling factors for atoms (example scaling factors for Te and W atoms)
% atom_factors = [0.00000044, 1.31473454, 6.35624338;
%               1.74897061, 4.43335644, 5.73543738;
%              -0.00000002, 5.44760823, 9.92929379;
%               1.74897108, 2.22834683, 9.30776171;
%               1.74897050, 6.06754902, 7.93340128;
%               0.00000040, 3.81415493, 7.73232446];
% 
% unit_cell_length=[3.49794  6.33827 15.4319];
% ax =  3.49794;
% by =  6.33827;
% cz =  15.4319;
% 
% % % % % % Initialize atomic positions
% main_atomic_position = [];
% neighbor_atomic_position=[];
% % Generate positions for Te and W atoms
% for i = 1:length(main_unit)
% 
%        for j = 1:size(atom_factors, 1)
%         main_atomic_position(end+1, :) = main_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
%        neighbor_atomic_position(end+1, :) = neighbor_unit(i, :) .*unit_cell_length+ atom_factors(j, :);
%        end
% 
% end
% 
% 
% %% device atomic position 
% Device_atomic_position=[];
% for i = 1:length(Device_unit)
% 
%        for j = 1:size(atom_factors, 1)
%          Device_atomic_position(end+1, :) = Device_unit(i, :) .*unit_cell_length+ atom_factors(j, :);  % Scale main unit for Te atoms
% 
%        end
% 
% end
% 
% %%%% next calculation 
% 
% Device_atomic_position=Device_atomic_position;
% 
% finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
% [number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);
% 
% x_uniq=unique(finalc_w_extra_lrp(:,1));
% y_uniq=unique(finalc_w_extra_lrp(:,2));
% z_uniq=unique(finalc_w_extra_lrp(:,3));
% xyzlengthuh=[(max(x_uniq)-min(x_uniq)) (max(y_uniq)-min(y_uniq)) (max(z_uniq)-min(z_uniq))]
% Vds=0.2;%%%drain voltage
% Vgs=[0 0.3 0.6 1.2];%%%%%gate voltage
% 
% Unit_cell_bias_drain=zeros(44);
% Unit_cell_bias_source=zeros(44);
% for i=1:44
%     Unit_cell_bias_drain(i,i)=Vds;
%     Unit_cell_bias_source(i,i)=0;
% 
% end
% 
% BBBB=diag(ones(1,num_cell));
% big_H_bias_drain=kron(BBBB,Unit_cell_bias_drain);
% big_H_bias_source=kron(BBBB,Unit_cell_bias_source);
% 
% 
% Length_Device=max(x_uniq)-min(x_uniq);
% potential_across_length=zeros(44*num_cell*Unit_Cells,44*num_cell*Unit_Cells);
% potential_across_hight=zeros(number_atoms,1);
% oooo=0;
% for iii=1:number_atoms
%        hi=Vds*(finalc_w_extra_lrp(iii,1)-min(x_uniq))/Length_Device;
%        if(mod(iii,6)==1||mod(iii,6)==4)
%        for jjj=1:6
%            potential_across_length(jjj+oooo,jjj+oooo)=hi;
% 
%        end
%        oooo=oooo+6;
%        end
% 
%        if(mod(iii,6)==2||mod(iii,6)==3)
%        for jjj=1:6
%            potential_across_length(jjj+oooo,jjj+oooo)=0;
% 
%        end
%        oooo=oooo+6;
%        end
% 
% 
% 
% 
%        if(mod(iii,6)==5||mod(iii,6)==0)
%        for jjj=1:10
%            potential_across_length(jjj+oooo,jjj+oooo)=hi;
%        end
%        oooo=oooo+10;
%        end
% end
% 
% 
% % % mod(8,6)
% % % % [xxGrid, yyGrid]=meshgrid(x_uniq ,y_uniq);
% 
% 
% % % % % % \\\\\\\\\\\\\\\\\ Atom position //////////////////
% na = num_cell/3;
% at1=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
% at=[1 2 3 4 5 6;1 2 3 4 5 6;1 2 3 4 5 6];
% for j=1:na-1
%     at=[at;at1];
% 
% end
% 
% % \\\\\\\\\\\\\\\ Calling the function///////////////
% cell_pos = cp;
% cell_pos1 = cpp;
% atomst = at;
% atomst1 = at;
% proj = [6 6 6 6 10 10];
% HNN = 5;
% Vext = 0;
% [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
% 
% 
% % ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
% %%% Device calculation
% eta=1e-8;
% count1=44*num_cell;
% count2=count1*Unit_Cells; 
% tox = 7; %in angstorm(oxide thickness)
% eps_ox = 25; % relative permittivity of HfO2
% 
% 
% 
% [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
% H= hamiltonian_TMD_device_H(bigH,bigU,Unit_Cells);
% 
% % /////// device calculation \\\\\\\\\\\
% Ef = 3.13195; % Fermi energy of Bulk Bi(111) plane 
% kT = 0.0259;% at 300 K
% eta = 1e-6;
% tol = 0.001;
% 
% 
% h = 6.626e-34;
% 
% q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h);
% I0=(2*q^2)/h;
% I0=(q^2)/hbar/2/pi;
% 
dTe1  =       6.35624338;
dTe2   =      5.73543738;
dTe3    =     9.92929379;
dTe4     =    9.30776171;
dW1      =   7.93340128;
dW2        =  7.73232446;
% % % Temperature
% Temp=300;
% count=0
% % % Define Fermi function
% f=@(E,Temp)1./(1+exp(E/(kb*Temp)));
% 
% 
% 
% 
% 
% 
% %% ///////////////// Starting Calculation \\\\\\\\\\\\\\\\\
% figure(1)
% 
% for kk=1:length(Vgs)
% 
% %%%%%% External voltage apply 
% [bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);%%%%% [hr,hc]=size(bigH); 
% Voltage=Vgs(kk);
% 
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
% 
BBBB=diag(ones(1,num_cell*Unit_Cells));
big_H_gate_voltage=kron(BBBB,Unit_cell_H);
big_H_exeternal_voltage=big_H_gate_voltage;
% UL=potential_across_length;
% U=UL;
% Phi=diag(UL);
% 
% 
%     mu1=Ef;
%     mu2=Ef-Vds;
%     E = linspace(mu2,mu1,20);
%     dE=E(3)-E(2);
%     delphi = 1;
%     Ite = 0;
%     f1=1./(1+exp((E-mu1)./kT));
% 	f2=1./(1+exp((E-mu2)./kT));
%     epsilon=1;
%     rho = 0;
% while (delphi>tol)
%     rho = 0;
%     Phi_old = Phi;
% 
% for k1=1:length(E)
% h=bigH;
% ud=bigUbar;
% u=bigU;
% H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
% [G1,countG1]=lead_eps(h,ud,eta,E(k1),count1,count2);
% [G2,countG2]=lead_eps(h+big_H_bias_drain,u,eta,E(k1),count1,count2);
% sig1=ud*G1*u;
% sig2=u*G2*ud;
% 
% sig1f=zeros(count2);
% sig2f=sig1f;
% for m=1:count1
%     for n=1:count1
%         sig1f(m,n)=sig1(m,n);
%         sig2f(m+count2-count1,n+count2-count1)=sig2(m,n);
%     end
% end
% 
% gam1=1i*(sig1f-sig1f');
% gam2=1i*(sig2f-sig2f');
% sigmain=gam1*f1(k1)+gam2*f2(k1);
% sigmaout=gam1*(1-f1(k1))+gam2*(1-f2(k1));
% 
% Gr=inv((E(k1)*eye(count2))-H-big_H_exeternal_voltage-Phi-sig1f-sig2f);
% Ga=Gr';
% Gn=Gr*sigmain*Ga;
% Gp=Ga*sigmaout*Gr;
% rho=rho+(1/(2*pi))*Gn*dE/(ax*1e-10*by*1e-10);
% T(k1)=trace(gam1*Gr*gam2*Ga);
% 
% 
% end
%   prho = real(diag(rho));
% 
%   %%% new prho 
%   %  every wte2 has [6 6 6 6 10 10] orbital
% jjjj=1;
% for iiii=1:44:length(prho)
%     prho_new(jjjj)=sum(prho(iiii:iiii+6-1));
%     prho_new(jjjj+1)=sum(prho(iiii+6:iiii+12-1));
%     prho_new(jjjj+2)=sum(prho(iiii+12:iiii+18-1));
%     prho_new(jjjj+3)=sum(prho(iiii+18:iiii+24-1));
%     prho_new(jjjj+4)=sum(prho(iiii+24:iiii+34-1));
%     prho_new(jjjj+5)=sum(prho(iiii+34:iiii+44-1));
%     jjjj=jjjj+6;
% end
% 
% 
% 
% 
% 
%   BCconduh={'Dirichlet','Dirichlet'};
% 
%     BCVOuh=[0 Vds Voltage -Voltage];
%     % function [UUUet,UUet,modeletn,uet,f111Det,resultset]=MATLAB_PDE_3D(finalcuh,f1Duh,BCVOuh,BCconduh)
%     [UUUet,fcord,UUet,modeletn,uet,f111Det,resultset,potential_new_input]=PoissonJS3D(finalc_w_extra_lrp,prho_new,BCVOuh,BCconduh,count1,count2,Vds,xyzlengthuh);
%   Phi=potential_new_input;
%   delphi = max(abs(Phi-Phi_old))
% %   Phi1=potential_new_input;
%  Phi = 0.25.*(Phi)+0.75.*Phi_old;
% % %     Phi = Phi1;
%     % delphi = max(abs(Phi-Phi_old))
%     Ite = Ite+1
%     if Ite>50
%         break;
%     end
% 
% end
% 
% dI=IE*(f(E-Ef,Temp)-f(E-Ef+Vds,Temp)).*T*dE;
% I(kk)=sum(dI);
%     count=count+1
% potential=UUet;
% 
% 
% 
% % % save Igsvsvgs I
% % % hold on
% % % h=plot(Vgs,I*10^6,'b');
% % % set(h,'linewidth',[2.0])
% % % set(gca,'Fontsize',[24])
% % % grid on
% % % xlabel(' Voltage Vgs ( V ) --> ')
% % % ylabel(' Current Ids (µA) --> ')
% 
% % 
% % % % % % % Extract and scale coordinates
% x = finalc_w_extra_lrp(:, 1); % x-coordinates
% y = finalc_w_extra_lrp(:, 2); % y-coordinates
% z = finalc_w_extra_lrp(:, 3); % z-coordinates
% Ec = -UUet; % Potential values (negative of UUet)
% 
% 
% % Fixed z-coordinate
% z_fixed = 9.30776171;
% 
% % Ensure unique values for x and y (grid points)
% x_uniq = unique(x);
% y_uniq = unique(y);
% 
% 
% 
% % Interpolate potential values at given coordinates
% % % V_filtered = griddata(x, y, z, UUet,coordinate_final(:, 1), coordinate_final(:, 2), coordinate_final(:, 3), 'nearest');
% 
% % Filter for fixed z-coordinate
% mask = (z == z_fixed); % Tolerance for fixed z
% x_filtered = x(mask); % Filtered x-coordinates
% y_filtered = y(mask); % Filtered y-coordinates
% V_filtered_z = Ec(mask); % Filtered potential values at z_fixed
% for ii=1:num_cell
%     x_last_1(ii,1)=max(x_uniq);
%     y_last_1(ii,1)=y_filtered(ii,1);
%     V_last_1(ii,1)=-Vds;
% end
% for ii=1:num_cell
%     x_last_2(ii,1)=max(x_uniq)+(x_filtered(4)-x_filtered(3));
%     y_last_2(ii,1)=y_filtered(ii,1);
%     V_last_2(ii,1)=-Vds;
% end
% 
% 
% for ii=1:num_cell
%     x_first_1(ii,1)=min(x_uniq)-(x_filtered(4)-x_filtered(3));
%     y_first_1(ii,1)=y_filtered(ii,1);
%     V_first_1(ii,1)=0;
% end
% % 
% % 
% % % x_filtered=[x_first_1; x_filtered; x_last_1;x_last_2];
% % % y_filtered=[y_first_1; y_filtered; y_last_1;y_last_2];
% % % V_filtered_z=[V_first_1; V_filtered_z; V_last_1;V_last_2];
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %%% gn plot of 2nd band 
% % % [fitresult, gof] = createFit(x_filtered(:,1),y_filtered(:,1),V_filtered_z);
% % % axis tight
% % % zlabel('Potential');
% % % movegui('west');
% % % xlabel('width');
% % % ylabel('length');
% % % % 
% % % % 
% % 
% % 
% % 
% % %%%%% 2 D plot 
% 
% 
% x_filtered=[x_first_1; x_filtered; x_last_1];
% y_filtered=[y_first_1; y_filtered; y_last_1];
% V_filtered_z=[V_first_1; V_filtered_z; V_last_1];
% 
% X_unique=unique(x_filtered);
% pp=1;
% V_unique=diag(zeros(length(X_unique)));
% for ii=1:length(X_unique)
%     V_unique(ii,1)=V_filtered_z(pp,1);
%     pp=pp+num_cell;
% end
% 
% plot(X_unique, V_unique, 'b', 'LineWidth', 1.9); % Plot with blue color and specified line width
% xlabel('X direction'); % Label the x-axis (modify as needed)
% ylabel('Potential'); % Label the y-axis (modify as needed)]
% title('Energy band diagram'); % Title of the plot
% grid on; % Turn on the grid for better visualization
% xlim([min(X_unique), max(X_unique)]); % Tighten only the x-axis
% hold on
% end
% 
% % % save Igsvsvgs I
% % % hold on
figure(1)
h=plot(Vgs,real(I)*10^6,'b');
set(h,'linewidth',[2.0])
set(gca,'Fontsize',[24])
grid on
xlabel(' Voltage Vgs ( V ) --> ')
ylabel(' Current Ids (µA) --> ')