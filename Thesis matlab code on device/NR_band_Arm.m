clc
clear 
close all
caption = 'wte2';
%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR =8;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........
Ef=3.13195;
% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
xa(1)=0;xa(2)=1;xa(3)=2;
ya(1)=0;ya(2)=0;ya(3)=0;
for j=0:num_cell-4
    
    xa(4+j)=xa(1+j)+3;
    ya(4+j)=ya(1+j)+0;
end
za=zeros(1,length(xa))';
cp = [xa',ya',za];

% /////////////// Neighbour Unit \\\\\\\\\\\\\\\\\\\\\\\\
xb(1)=0;xb(2)=1;xb(3)=2;
yb(1)=1;yb(2)=1;yb(3)=1;
for j=0:num_cell-4
    
    xb(4+j)=xb(1+j)+3;
    yb(4+j)=yb(1+j)+0;
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
cell_pos = cp;%%%% unit cell
cell_pos1 = cpp;%%% neighbour unit cell
atomst = at;%%% unit cell atoms
atomst1 = at;%%% neighbour cell atoms
proj = [6 6 6 6 10 10];%%%%
HNN = 3;
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
dTe1  =       6.35624338;
dTe2   =      5.73543738;
dTe3    =     9.92929379;
dTe4     =    9.30776171;
dW1      =   7.93340128;
dW2        =  7.73232446;

Vgs =[0]; % Top or Bottom Gate Voltage;


Vds=0.0;

 Unit_cell_bias=zeros(44);
for i=1:44
    Unit_cell_bias(i,i)=Vds;
    
end


%%%%%gate voltage
ts=3.3;eps_s=1.6;ti=3.17;eps_i=3.29;width=4.59; kb=8.3333e-05;Temp=300;% % Temperature
kT = 0.0259;
h_plank = 6.626e-34;
q = 1.6e-19; hbar = 1.06e-34; zplus=1i*1e-3; kb=8.3333e-05; IE=(2*q*q/ h_plank);





for ii=1:length(Vgs)

%%%%%% External voltage apply 
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
Voltage=Vgs(ii);
  Electric_field=-(2*Voltage)/(dTe2-dTe3);
   Unit_cell_H=zeros(44);



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
big_H_bias_voltage=kron(BBBB,Unit_cell_bias);


bigH=bigH+big_H_exeternal_voltage+big_H_bias_voltage;


%% External voltage apply end 
        
        HB=cell(3,1);
        HB{1,1}=sparse(bigUbar);
        HB{2,1}=sparse(bigH);
        HB{3,1}=sparse(bigU);
        unique_trans_vec=[0 -1 0;0 0 0;0 1 0];
        [hr,hc]=size(HB{2,1});
        diag(HB{2,1});
        
        Es=[];
        [kpr,kpc]=size(kpoints);
        xx=1
        yy=1
        for k=1:kpr
            Hk=zeros(hr,hc);
            for j=1:3
                Hb=HB{j,1};
                Hk=Hk+HB{j,1}.*exp(1i*sum(unique_trans_vec(j,:).*kpoints(k,:)));
            end
            [W,Ek]=eig(Hk);
                 En=sort(real(diag(Ek)));
                Es(:,k)=En-Ef;
                band_diagram_iteration_left = kpr-k
            
            
        end
        
        Ess=Es';

        [row,column]=size(Ess)

        for i=1:column
            if(mod(i,2)==1)
                Esup(:,xx)=Ess(:,i);
                xx=xx+1;
            end
            if(mod(i,2)==0)
                Esdown(:,yy)=Ess(:,i);
                yy=yy+1;
            end

        end


       
figure
plot(Ess,'b','LineWidth',1.0)
hold on
% plot(Ess(:,587),'r','LineWidth',1.2)
% plot(Ess(:,588),'r','LineWidth',1.2)
% plot(Ess(:,589),'r','LineWidth',1.2)
% plot(Ess(:,590),'r','LineWidth',1.2)
axis tight
 ylim([-0.4 +0.4])
 hold on
plot(0*ones(kpr),'k--','LineWidth',0.9)
hold off
% Set the y-axis and x-axis labels
ylabel('E-EF (eV)');
xlabel('Wave Vector');
end
