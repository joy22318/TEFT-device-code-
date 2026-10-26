function [delCube, Dz, Lz] = Diffmat_1D_Poisson(z, nz, addfe)
if nargin < 3
    addfe = 0;
end

delz = abs(gradient(z));

% First and second derivative operators for the z-direction
Dz = spdiags([-ones(nz,1), 0*ones(nz,1), ones(nz,1)], [-1, 0, 1], nz, nz); % First derivative
D2z = spdiags([-ones(nz,1), 2*ones(nz,1), -ones(nz,1)], [-1, 0, 1], nz, nz); % Second derivative

if addfe % Taking outer boundary points into account
    Dz(1,1:3) = [-3 4 -1];
    Dz(end,end-2:end) = [1 -4 3];
    
    D2z(1,1:4) = [2 -5 4 -1];
    D2z(end,end-3:end) = [-1 4 -5 2];
end

% Normalize the derivative matrices by the grid spacing
for j = 1:nz
   Dz(j,:) = Dz(j,:) * 1/(2 * delz(j)); % Final matrix operator for first derivatives 
   D2z(j,:) = D2z(j,:) * 1/(delz(j) * delz(j)); % Final matrix operator for second derivatives 
end

% Volume element in 1D (here it's just the spacing in z-direction)
delCube = delz(:);

% Laplacian operator in z-direction (since we're only focusing on z)
Lz = D2z;

end
