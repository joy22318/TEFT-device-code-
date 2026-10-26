function[kpoints_cart,kpoints_cryst]=kpgen(kpath,nopot)
kpathp=length(kpath(:,1));
kpoints_cryst=[];
for k=1:3
    kxf=[];
    for j=1:kpathp-1
        if j==1
            if kpath(j,k)==kpath(j+1,k)
                kx=kpath(j,k)*ones(1,nopot(j));
            else
                kx=linspace(kpath(j,k),kpath(j+1,k),nopot(j));
            end
            kxf=vertcat(kxf,kx');
        else
            if kpath(j,k)==kpath(j+1,k)
                kx=kpath(j,k)*ones(1,nopot(j));
            else
                kx=linspace(kpath(j,k),kpath(j+1,k),nopot(j));
            end
            kxn=kx(2:end);
            kxf=vertcat(kxf, kxn');
        end
    end
    kpoints_cryst(:,k)=kxf;
end
kpoints_cart=2*pi.*kpoints_cryst;
end
