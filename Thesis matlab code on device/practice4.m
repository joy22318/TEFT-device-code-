clc
clear 
close all
caption = 'wte2';
% ////////Creating the input\\\\\\\\\\\\\\\\\\\\\\
ZBNR = 4;  % Zigzag Bi nano ribbon , possible values = 4,8,12,16,20,24,28,........
% % % Calculating number of cells required
nn = ((ZBNR-4)/4)+1;
num_cell = ZBNR-nn; %number of cell possible values = 3,6,9,12,15,18,21,25,........

% % % //////// Main Unit \\\\\\\\\\\\\\\\\\\\\\
Unit_Cells=5;
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

x_uniq=unique(finalc_w_extra_lrp(:,1));
y_uniq=unique(finalc_w_extra_lrp(:,2));
z_uniq=unique(finalc_w_extra_lrp(:,3));
xyzlengthuh=[(max(x_uniq)-min(x_uniq)) (max(y_uniq)-min(y_uniq)) (max(z_uniq)-min(z_uniq))]
Vd=0.1;%%%drain voltage
Vs=-0.1;
Vds=0.2;
Vgs=0.3;%%%%%gate voltage


Unit_cell_bias_drain=zeros(44);
Unit_cell_bias_source=zeros(44);
for i=1:44
    Unit_cell_bias_drain(i,i)=Vd;
    Unit_cell_bias_source(i,i)=Vs;

end
BBBB=diag(ones(1,num_cell));
big_H_bias_drain=kron(BBBB,Unit_cell_bias_drain);
big_H_bias_source=kron(BBBB,Unit_cell_bias_source);


Length_Device=max(x_uniq)-min(x_uniq);
potential_across_length=zeros(6*num_cell*Unit_Cells,1);
potential_across_hight=zeros(number_atoms,1);
oooo=0;
for iii=1:length(potential_across_length)
       hi=Vds*(finalc_w_extra_lrp(iii,1)-min(x_uniq))/Length_Device+Vs;
       
      
           potential_across_length(iii,1)=hi;

       
      
end

oooo=1;
for iii=1:number_atoms
       if(mod(iii,6)==1||mod(iii,6)==4)

           potential_across_height(oooo,1)=0;


       oooo=oooo+1;
       end 

       if(mod(iii,6)==2)

           potential_across_height(oooo,1)=+Vgs;


       oooo=oooo+1;
       end 


       if(mod(iii,6)==3)

           potential_across_height(oooo,1)=-Vgs;


       oooo=oooo+1;
       end 

       if(mod(iii,6)==5||mod(iii,6)==0)

           potential_across_height(oooo,1)=0;

       oooo=oooo+1;
       end


end



% Extract and scale coordinates
x = finalc_w_extra_lrp(:, 1); % x-coordinates
y = finalc_w_extra_lrp(:, 2); % y-coordinates
z = finalc_w_extra_lrp(:, 3); % z-coordinates
Ec = potential_across_length; % Potential values (negative of UUet)

% Fixed z-coordinate
z_fixed = 7.73232446;

% Ensure unique values for x and y (grid points)
x_uniq = unique(x);
y_uniq = unique(y);

% Interpolate potential values at given coordinates
% % V_filtered = griddata(x, y, z, UUet,coordinate_final(:, 1), coordinate_final(:, 2), coordinate_final(:, 3), 'nearest');

% Filter for fixed z-coordinate
mask = (z == z_fixed); % Tolerance for fixed z
x_filtered = x(mask); % Filtered x-coordinates
y_filtered = y(mask); % Filtered y-coordinates
V_filtered_z = Ec(mask); % Filtered potential values at z_fixed

% Visualize filtered potential
figure;
scatter3(x_filtered, y_filtered, V_filtered_z, 40, V_filtered_z, 'filled');
xlabel('Width (x)');
ylabel('Length (y)');
zlabel('Potential (V)');
title('Potential Distribution at Fixed z-coordinate');
colorbar;
axis tight;

% Optional: Fit the potential surface for visualization
[fitresult, gof] = createFit(x_filtered, y_filtered, V_filtered_z);

% Plot the fitted surface
figure;
plot(fitresult);
xlabel('Width (x)');
ylabel('Length (y)');
zlabel('Potential (V)');
title('Fitted Potential Surface');
colorbar;
axis tight;







%%%% adding points 


% Define constants and inputs
z_fixed = 7.73232446; % Fixed z-coordinate
z_new = z_fixed; % Optional: Keep z constant or adjust if needed

% Generate finer grid for x and y
x_range = linspace(min(x_filtered), max(x_filtered), 100); % Increase points in x-direction
y_range = linspace(min(y_filtered), max(y_filtered), 100); % Increase points in y-direction
[X_grid, Y_grid] = meshgrid(x_range, y_range); % Create a grid of x and y

% Interpolate the potential values for the finer grid at fixed z
V_interpolated = griddata(x_filtered, y_filtered, V_filtered_z, X_grid, Y_grid, 'cubic');

% Visualize the interpolated potential values
figure;
surf(X_grid, Y_grid, V_interpolated, 'EdgeColor', 'none'); % Smooth surface plot
xlabel('Width (x)');
ylabel('Length (y)');
zlabel('Potential (V)');
title('Interpolated Potential Distribution');
colorbar;
axis tight;

% Optional: Add the constant z value (if visualizing in 3D space)
Z_grid = z_new * ones(size(X_grid)); % Create constant z-plane
figure;
scatter3(X_grid(:), Y_grid(:), Z_grid(:), 20, V_interpolated(:), 'filled');
xlabel('Width (x)');
ylabel('Length (y)');
zlabel('Potential (V)');
title('Interpolated Potential Values at Fixed z-coordinate');
colorbar;
axis tight;

% Optional: Surface fitting with finer grid
[fitresult, gof] = createFit(X_grid(:), Y_grid(:), V_interpolated(:)); % Fit with finer data

% Plot fitted surface
figure;
plot(fitresult);
xlabel('Width (x)');
ylabel('Length (y)');
zlabel('Potential (V)');
title('Fitted Potential Surface on Finer Grid');
colorbar;
axis tight;







