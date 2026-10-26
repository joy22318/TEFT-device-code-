function [potential_new_input,finalc,potential]=poison_JS_horizon(ZBNR,Unit_Cells,prho,vds)

%% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = ZBNR;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........
% //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
%% ////////Declaration of Variable\\\\\\\\\\\\\\\\\\\\\\
%%%% Device calculation
% Unit_Cells = Unit_Cells; %Number of Unit Cells
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


x_max=max(Device_atomic_position(:,1))
x_min=min(Device_atomic_position(:,1))


%%%%% next calculation 

Device_atomic_position=Device_atomic_position.*10^-10;

finalc_w_extra_lrp = [Device_atomic_position(:,1),Device_atomic_position(:,2),Device_atomic_position(:,3)];
[number_atoms co_ordinate_numb]=size(finalc_w_extra_lrp);

%% ////////// Adding extra points \\\\\\\\\\\\\\\\\\\\\\\\
 k11=1;
 % % % ch_r = linspace(x_min,x_max,40);

%%%  % ch_r = 5.73543738:dz:9.92929379;
 ch_r=linspace(x_min*10^-10,x_max*10^-10,40);
 dz=ch_r(2)-ch_r(1);
% 
 n_extra = length(ch_r);
for i=1:length(finalc_w_extra_lrp)

    for j=0:n_extra-1
    finalc(k11+j,:)=finalc_w_extra_lrp(i,:);
    end
    k11=k11+n_extra;
end

for i=1:n_extra:length(finalc)
    for j=1:n_extra
        finalc(i+j-1,1)=ch_r(j);
    end

end
% 
% 
% 
% 
% %% Calling the system differential matrix function 
% 
% 
%  every wte2 has [6 6 6 6 10 10] orbital
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
f1D = griddata(A_X1(:,1),A_Y1(:,1),A_Z1(:,1),prho_new',finalc(:,1),finalc(:,2),finalc(:,3),'nearest');

eps0 = 8.854e-12;
eps = 15*eps0;
f1D = (f1D.*1e-19)./(eps); % Check again



%% Calling Differential matrix 
[H_poison] = Diffmat_1D( n_extra, number_atoms,dz);

% 
% 
% 
for ii=1:length(f1D)
     if(finalc(ii,3)==x_max*10^-10)
         f1D(ii)=vds;
     end
     if(finalc(ii,3)==x_min*10^-10)
         f1D(ii)=0;
     end

end
% 
potential =inv(H_poison)*(f1D);
% 
for ss=1:(length(potential)/n_extra)
    distance_etter=100;
    queryX = A_X1(ss, 1);
    queryY = A_Y1(ss, 1);
    queryZ = A_Z1(ss, 1);

    start=(ss-1)*n_extra+1;

    for kk=start:start+n_extra-1
    distances_cal = sqrt((finalc(kk,1) - queryX)^2 + (finalc(kk,2) - queryY)^2 + (finalc(kk,3) - queryZ)^2);
    if(distances_cal<distance_etter)
    potential_new(ss)=potential(kk);
    end

    distance_etter=distances_cal;
    end
end

hhhh=0;

% 
for ll=1:6:length(potential_new)

      for gg=1:6
      potential_new_input(gg+hhhh,1)=potential_new(ll);
      end
      for gg=7:12
      potential_new_input(gg+hhhh,1)=potential_new(ll+1);
      end
      for gg=13:18
      potential_new_input(gg+hhhh,1)=potential_new(ll+2);
      end
      for gg=19:24
      potential_new_input(gg+hhhh,1)=potential_new(ll+3);
      end
      for gg=25:34
      potential_new_input(gg+hhhh,1)=potential_new(ll+4);
      end
      for gg=35:44
      potential_new_input(gg+hhhh,1)=potential_new(ll+5);
      end
      hhhh=hhhh+44;
end
% 
% 
% 

