function [H_poison] = Diffmat_1D( n_extra, number_atoms,dz)

N=n_extra;
 H=eye(N,N);
for i=2:N-1
    H(i,i-1)=1;
    H(i,i+1)=1;
    H(i,i)=-2;
end

H(N,N)=dz^2;
H(1,1)=dz^2;

H=(1/dz^2).*H;

I=diag(ones(1,number_atoms));

H_poison=kron(I,H);
