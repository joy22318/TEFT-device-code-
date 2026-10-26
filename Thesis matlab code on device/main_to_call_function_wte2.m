clc
clear 
close all
caption = 'wte2';
%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR =32;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........
Ef=3.4055;
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

for ii=1:length(Vgs)

%%%%%% External voltage apply 
[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH); 
Voltage=Vgs(ii);
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
                Es(:,k)=En;
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
plot(kpoints(:,2),Esup,'b','LineWidth',1.9)
hold on
plot(kpoints(:,2),Esdown,'r--','LineWidth',1.9)
axis tight
 ylim([2.5 5])
 hold on
plot(kpoints(:,2),Ef*ones(kpr),'k--')
hold off
end
%% transport calculation 
% Unit_Cells=10;
% E2=3.4055;
% eta=1e-8;
% count1=44*num_cell;
% count2=count1*Unit_Cells;   
% f1=1;
% f2=0;
% %[bigH,bigU,bigUbar]=cal_H_U_Ubar(caption,cell_pos,cell_pos1,atomst,atomst1,proj,HNN,Vext);% [hr,hc]=size(bigH);
% %H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
% h=bigH;
% ud=bigUbar;
% u=bigU;
% H= hamiltonian_TMD_device_H(h,u,Unit_Cells);
% [G1,countG1]=lead_eps(h,ud,eta,E2,count1,count2);
% [G2,countG2]=lead_eps(h,u,eta,E2,count1,count2);
% sig1=ud*G1*u;
% sig2=u*G2*ud;
% 
% % % Calculating and plotting Gn
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
% 
% ppp=2
% gam1=1i*(sig1f-sig1f');
% gam2=1i*(sig2f-sig2f');
% sigmain=gam1*f1+gam2*f2;
% sigmaout=gam1*(1-f1)+gam2*(1-f2);
% 
% Gr=inv((E2*eye(count2))-H-sig1f-sig2f);
% Ga=Gr';
% Gn=Gr*sigmain*Ga;
%     Gp=Ga*sigmaout*Gr;
% 
% 
%  T=trace(gam1*Gr*gam2*Ga);
% 
%  for j=1:44
% 
%  p=j;
%   for i=1:(count2/44)
%       gn1(i,i)=Gn(p,p);
%       p=p+44;
% 
%   end
% 
% figure 
% Gnplot2(gn1,count1/44,count2/44,E2)
%  end
% 
% 
