clc
clear 
close all
A = importdata('wte2.projbands.txt');
% A = importdata('Bi_soc1.txt')
%%
x=1;
for i=1:110:8910
    k(x)=A(i,2);

for j=1:110
   
    E(j,x)=A(i+j-1,3);
end
 

for j=1:110
    sum=0;
   for l=85:94
       sum=sum+A(i+j-1,l)+sum+A(i+j-1,l+20);
   
   end
   d(j,x)=sum;
end

for j=1:110
    sum=0;
   for l=17:22
       sum=sum+A(i+j-1,l)+A(i+j-1,l+18)+A(i+j-1,l+36)+A(i+j-1,l+54);
   
   end
   p(j,x)=sum;
end
    
    
x=x+1;
end
% figure()
% plot(k,Es1,k,Es2,k,Es3,k,Es4,k,Es5,k,Es6,k,Es7,k,Es8)
%% 
figure()
for i=1:110
patch([k nan],[E(i,:) nan],[d(i,:) nan],[d(i,:) nan],'edgecolor','interp');
colorbar;colormap(jet)
title('Projected Band structure on d orbital for W')
xlabel('Wave number k')
ylabel('Energy (eV)')
set(get(gca,'xlabel'),'FontWeight','Bold');
set(get(gca,'ylabel'),'FontWeight','Bold');
movegui('east');
ylim([-6,9])
hold on
end
hold off

figure()
for i=1:110
patch([k nan],[E(i,:) nan],[p(i,:) nan],[p(i,:) nan],'edgecolor','interp');
colorbar;colormap(jet)
title('Projected Band structure on P orbital for Te')
xlabel('Wave number k')
ylabel('Energy (eV)')
set(get(gca,'xlabel'),'FontWeight','Bold');
set(get(gca,'ylabel'),'FontWeight','Bold');
movegui('east');
  ylim([-6,9])
hold on
end
hold off


