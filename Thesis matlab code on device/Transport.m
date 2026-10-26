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




%% transport calculation 
Unit_Cells=16;
E2=Ef;
eta=1e-8;
count1=44*num_cell;
count2=count1*Unit_Cells;   
f1=0;
f2=1;
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


 T=trace(gam1*Gr*gam2*Ga);

 for j=1:44

 p=j;
  for i=1:(count2/44)
      gn1(i,i)=Gn(p,p);
      p=p+44;

  end

figure 
Gnplot2(gn1,count1/44,count2/44,E2)
axis tight
zlabel('\mid\Psi\mid^{2}');
movegui('west');
xlabel('Length --->')
ylabel('Width --->')
ax=gca;
ax.FontSize=10;
ylh = get(gca,'ylabel');
set(ylh, 'Rotation',-25,'VerticalAlignment','top', 'HorizontalAlignment','right')
xlh = get(gca,'xlabel');
set(xlh, 'Rotation',15,'VerticalAlignment','top', 'HorizontalAlignment','left')
set(get(gca,'zlabel'))
 end


