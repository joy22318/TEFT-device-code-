function [UUUet,fcord,UUet,modeletn,uet,f111Det,resultset,potential_new_input]= PoissonJS3D(finalcuh,f1Duh,BCVOuh,BCconduh,count1,count2,Vds,xyzlengthuh)
% % uXXX,uYYY,uZZZ,finalcordxxx,xyzlengthuh,
  ang=1e-10;
  ev = 1.6e-19; 
  epsilon0=8.85410e-12;
  epsilonch=25;

  uXXX=unique(finalcuh(:,1));
  uYYY=unique(finalcuh(:,2));
  uZZZ=unique(finalcuh(:,3));
        


        [xxGrid, yyGrid,zzGrid]=meshgrid(uXXX*ang, uYYY*ang, uZZZ*ang);
	    
        XT=xxGrid(:);
        YT=yyGrid(:);
        ZT=zzGrid(:);
        K=convhull(XT,YT,ZT);
        nodes=[XT YT ZT]';
        [~,finalcordxxx]=ismember(finalcuh.*ang,nodes','rows');
        elements=K';
        modeletn = createpde();
        geometryFromMesh(modeletn,nodes,elements);
        mesh=generateMesh(modeletn);
        incordd=modeletn.Mesh.Nodes;
        incordd=incordd';
        f1111D=griddata(finalcuh(:,1)*ang,finalcuh(:,2)*ang,finalcuh(:,3)*ang,f1Duh,incordd(:,1),incordd(:,2),incordd(:,3),'nearest');       

        f111Det=-f1111D.*ev/(epsilon0*epsilonch);
    
        
        % Voxx=(2*epsilonchuh*toxuh*BCVOuh(4)/(epsilonoxuh*xyzlengthuh(3)*angg))/(1+2*epsilonchuh*toxuh/(epsilonoxuh*xyzlengthuh(3)*angg));
        Vtopp=BCVOuh(4);%-Voxx;
        Vbottt=BCVOuh(3);%+Voxx;
        
        if strcmp(BCconduh{1},'Dirichlet')
            applyBoundaryCondition(modeletn,"dirichlet","Face",3,"u",BCVOuh(1));
            applyBoundaryCondition(modeletn,"dirichlet","Face",5,"u",BCVOuh(2));
        elseif strcmp(BCconduh{1},'Neumann')
            applyBoundaryCondition(modeletn,"neumann","Face",3,"g",BCVOuh(1)/(xyzlengthuh(2)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",5,"g",BCVOuh(2)/(xyzlengthuh(2)*ang));
        end
        if strcmp(BCconduh{2},'Dirichlet')
            applyBoundaryCondition(modeletn,"dirichlet","Face",6,"u",Vtopp);
            applyBoundaryCondition(modeletn,"dirichlet","Face",1,"u",Vbottt);  
        elseif strcmp(BCconduh{2},'Neumann')
            applyBoundaryCondition(modeletn,"neumann","Face",6,"g",Vtopp/(xyzlengthuh(3)*ang));
            applyBoundaryCondition(modeletn,"neumann","Face",1,"g",Vbottt/(xyzlengthuh(3)*ang));
        elseif strcmp(BCconduh{2},'Mixed')
            applyBoundaryCondition(modeletn,"dirichlet","Face",6,"u",BCVOuh(3));
            applyBoundaryCondition(modeletn,"neumann","Face",1,"g",BCVOuh(4)/(xyzlengthuh(3)*ang));
        elseif strcmp(BCconduh{2},'No_BC') 
        end
    
        f_fun = @(location,state)ffcoeffunction(location,incordd,f111Det); 
        specifyCoefficients(modeletn,"m",0,"d",0,"c",1,"a",0,"f",f_fun);
        generateMesh(modeletn);
        resultset = solvepde(modeletn);
        uet = resultset.NodalSolution;
        fcord=resultset.Mesh.Nodes;
        fcord=fcord';
      
        UUUet = griddata(fcord(:,1),fcord(:,2),fcord(:,3),uet,xxGrid(:),yyGrid(:),zzGrid(:),'linear');
        UUet=UUUet(finalcordxxx);%%%% very important output
        %%%% new potential profile 

hhhh=0;


for ll=1:6:length(UUet)
   
      for gg=1:6
      potential_new_input(gg+hhhh,1)=UUet(ll);
      end
      for gg=7:12
      potential_new_input(gg+hhhh,1)=UUet(ll+1);
      end
      for gg=13:18
      potential_new_input(gg+hhhh,1)=UUet(ll+2);
      end
      for gg=19:24
      potential_new_input(gg+hhhh,1)=UUet(ll+3);
      end
      for gg=25:34
      potential_new_input(gg+hhhh,1)=UUet(ll+4);
      end
      for gg=35:44
      potential_new_input(gg+hhhh,1)=UUet(ll+5);
      end
      hhhh=hhhh+44;
end



    

    	function fet = ffcoeffunction(locationet,finalcordet,f1Dett)
		    f1Dett=f1Dett';
		    temcord=[locationet.x' locationet.y' locationet.z'];
		    [tempr,~]=size(temcord);
            fet=zeros();
		    for i=1:tempr
			    % tempcr=temcord(i,3);
			    distances = pdist2(temcord(i,:),finalcordet);
			    % Find the minimum distance and its index for each point in coordinates_set1
			    [~, nearest_indices] = min(distances, [], 2);
			    fet(1,i)=f1Dett(nearest_indices);
		    end
        end

end





