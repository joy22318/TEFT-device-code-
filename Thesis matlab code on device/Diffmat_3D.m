function [delCube, Dx3D, Dy3D, Dz3D,L3D] = Diffmat_3D(x,y,z,nx,ny,nz,addfe)
if nargin<7
    addfe=0;
end

delx=abs(gradient(x));
dely=abs(gradient(y));
delz=abs(gradient(z));

[delX, delY,delZ]=meshgrid(delx,dely,delz);
delCubic=delX.*delY;
delCube=delCubic(:);

Dx = spdiags([-ones(nx,1), 0*ones(nx,1), ones(nx,1)], [-1, 0, 1], nx, nx);% first derivative
Dy = spdiags([-ones(ny,1), 0*ones(ny,1), ones(ny,1)], [-1, 0, 1], ny, ny);
Dz = spdiags([-ones(nz,1), 0*ones(nz,1), ones(nz,1)], [-1, 0, 1], nz, nz);

D2x = spdiags([-ones(nx,1), 2*ones(nx,1), -ones(nx,1)], [-1, 0, 1], nx, nx);% 2nd derivative
D2y = spdiags([-ones(ny,1), 2*ones(ny,1), -ones(ny,1)], [-1, 0, 1], ny, ny);
D2z = spdiags([-ones(nz,1), 2*ones(nz,1), -ones(nz,1)], [-1, 0, 1], nz, nz);

if addfe %%%%taking outer boundary points into account 
    Dx(1,1:3)=[-3     4    -1];
    Dx(end,end-2:end)=[1    -4     3];
    Dy(1,1:3)=[-3     4    -1];
    Dy(end,end-2:end)=[1    -4     3];
    Dz(1,1:3)=[-3     4    -1];
    Dz(end,end-2:end)=[1    -4     3];

    D2x(1,1:4)=[2 -5 4 -1];
    D2x(end,end-3:end)= [-1     4    -5     2];
    D2y(1,1:4)=[2 -5 4 -1];
    D2y(end,end-3:end)= [-1     4    -5     2];
    D2z(1,1:4)=[2 -5 4 -1];
    D2z(end,end-3:end)= [-1     4    -5     2];
end

for j=1:nx
   Dx(j,:)=Dx(j,:).*1/(2*delx(j));%%%%final matrix operator for first derivatives 
   D2x(j,:)=D2x(j,:).*1/(delx(j)*delx(j));%%%%final matrix operator for second derivatives 
end

for j=1:nz
   Dz(j,:)=Dz(j,:).*1/(2*delz(j));%%%%final matrix operator for first derivatives 
   D2z(j,:)=D2z(j,:).*1/(delz(j)*delz(j));%%%%final matrix operator for second derivatives 
end

for j=1:ny
   Dy(j,:)=Dy(j,:).*1/(2*dely(j));%%%%final matrix operator for first derivatives 
   D2y(j,:)=D2y(j,:).*1/(dely(j)*dely(j));%%%%final matrix operator for second derivatives 
end

Ix = speye(nx);
Iy = speye(ny);
Iz = speye(nz);   

Dx3D= kron(kron(Dx, Iy),Iz);
Dy3D= kron(Ix,kron(Dy, Iz));
Dz3D= kron(Ix,kron(Iy, Dz));

L3D= kron(Ix,kron(Iy, D2z)) + kron(Ix,kron(D2y, Iz)) + kron(kron(D2x, Iy),Iz);

end