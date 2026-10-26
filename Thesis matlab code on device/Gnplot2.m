%%% this function plot out zig zag nano ribbon 

function []=Gnplot2(Gn,endcount1,endcount2,energy)
psi2=diag(Gn);
sum(psi2)
for count1=1:endcount2
    res1=ceil(count1/endcount1);
    res2=mod(count1,endcount1);
    if res2==0
        res2=endcount1;
    end
    psi5(res2,res1)=psi2(count1);
end
%figure;
%str = sprintf('Energy = %f',energy);
surf(abs(psi5),'FaceColor','interp',...
'EdgeColor','none',...
'FaceLighting','phong')
%title(str);
%camlight left;