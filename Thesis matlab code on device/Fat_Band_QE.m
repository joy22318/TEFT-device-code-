function [numBut,xstut,ystut,wstut,ppput,ppputt]=Fat_Band_QE(input_file)
    config = read_input_file(input_file);

    write_prefix=config.WRITE_PDOS.prefix;
    write_inputdir=config.WRITE_PDOS.inputdir;
    Nspin_Type=config.WRITE_PDOS.Nspin_Type;
    write_outdir=config.WRITE_PDOS.outdir;
    Calprefix=config.WRITE_PDOS.Calprefix;
    E_Fermi=config.WRITE_PDOS.Efermi;
    Overwrite=config.WRITE_PDOS.Overwrite;
    Write_Delimiter=config.WRITE_PDOS.Delimiter;

    if ~isempty(write_inputdir)
        new_filename_band_in=sprintf('%s.%s.in',write_prefix,Write_Delimiter);
        capname_band_in=sprintf('%s/%s',write_inputdir,new_filename_band_in);
    end

    Want_States=config.PLOT_PDOS.Want_States;  % Cell array default
    Band_Limit=config.PLOT_PDOS.Band_Limit;
    Color_Matrix=config.PLOT_PDOS.Color_Matrix;  % Default matrix
    Sorting_Type=config.PLOT_PDOS.Sorting_Type;  % Cell array default
    xy_limit=config.PLOT_PDOS.xy_limit;
    Interpoint=config.PLOT_PDOS.Interpoint;  % Default matrix
    Weight_Multiplier=config.PLOT_PDOS.Weight_Multiplier;  % Cell array default
    Justify_Zero=config.PLOT_PDOS.Justify_Zero;
    Figure_Spec=config.PLOT_PDOS.Figure_Spec;  % Default matrix
    Draw_Fermi=config.PLOT_PDOS.Draw_Fermi;  % Cell array default
    Legend_Details=config.PLOT_PDOS.Legend_Details;
    Aspect_Ratio=config.PLOT_PDOS.Aspect_Ratio;  % Default matrix
    Export_Figure=config.PLOT_PDOS.Export_Figure;  % Default matrix
    Plot_Delimiter=config.PLOT_PDOS.Delimiter;
    % Figure_Caption=config.PLOT_PDOS.Figure_Caption;


    if isempty(write_prefix) && isempty(write_inputdir) && isempty(write_outdir) && isempty(Nspin_Type) && isempty(Calprefix) && isempty(E_Fermi) && isempty(Overwrite) 
        numBut=0;
        xstut=0;
        ystut=0;
        wstut=0;
    else
        [numBut,xstut,ystut,wstut,capnewut]=write_kp_band_psi(write_prefix,write_inputdir,Nspin_Type,write_outdir,Calprefix,E_Fermi,Overwrite);
    end

    if ~isfield(config.PLOT_PDOS,"prefix")
        plot_prefix=write_prefix;
    else
        plot_prefix=config.PLOT_PDOS.prefix;
    end
    if ~isfield(config.PLOT_PDOS,"inputdir")
        plot_inputdir=capnewut;
    else
        plot_inputdir=config.PLOT_PDOS.inputdir;
    end
    if ~isfield(config.PLOT_PDOS,"outdir")
        plot_outdir=capnewut;
    else
        plot_outdir=config.PLOT_PDOS.outdir;
    end

    
    if isfield(config.PLOT_PDOS,"inputdir")
        new_filename_band_in=sprintf('%s.%s.in',plot_prefix,Plot_Delimiter);
        capname_band_in=sprintf('%s/%s',plot_inputdir,new_filename_band_in);
    end

    if isempty(Want_States) && isempty(Band_Limit) && isempty(Color_Matrix) && isempty(Sorting_Type) && isempty(xy_limit) && isempty(Interpoint) && isempty(Weight_Multiplier) && isempty(Justify_Zero) && isempty(Figure_Spec) && isempty(Draw_Fermi) && isempty(Aspect_Ratio) && isempty(Export_Figure) && isempty(Legend_Details)
        ppput=0;
        ppputt=0;
    else
        [ppputt]=fat_band_plot_col_var(plot_prefix,plot_inputdir,Want_States,plot_outdir,Band_Limit,Color_Matrix,Sorting_Type,xy_limit,Interpoint,Weight_Multiplier,Justify_Zero,Figure_Spec,Draw_Fermi,Legend_Details,Aspect_Ratio,Export_Figure);
        [ppput]=fat_band_plot(plot_prefix,plot_inputdir,Want_States,plot_outdir,Band_Limit,Color_Matrix,Sorting_Type,xy_limit,Interpoint,Weight_Multiplier,Justify_Zero,Figure_Spec,Draw_Fermi,Legend_Details,Aspect_Ratio,Export_Figure);
    end


    function [numB,xstt,ystt,wstt,capnew]=write_kp_band_psi(varargin)
        newvararginnew = varargin;
        varargin = cell(1, 7); % Initialize with 7 empty cells
    
        for j = 1:length(newvararginnew)
            if isempty(newvararginnew{j})
                % If input is an empty array, store empty array
                varargin{j} = [];
            else
                % If input is zero or any other non-empty value, store it as is
                varargin{j} = newvararginnew{j};
            end
        end
        
        if nargin <8
            if isempty(varargin{1})
                error('Uable to read pdos file prefix.pdos.out');
            else
                prefix=varargin{1};
            end
            if isempty(varargin{2})
                error('Could find pdos file without input directory');
            else
                inputdir=varargin{2};
            end
            if isempty(varargin{3})
                error('Define nspin_type="no_soc" or "soc" or "spin_up_down"');
            else
                nspin_type=varargin{3};
            end
            if isempty(varargin{4})
                outdir='.';
            else
                outdir=varargin{4};
            end
            if isempty(varargin{5})
                calprefix=prefix;
            else
                calprefix=varargin{5};
            end
            if isempty(varargin{6})
                Efermi=1;
            else
                Efermi=varargin{6};
            end
            if isempty(varargin{7})
                overwrite=0;
            else
                overwrite=varargin{7};
            end
            
        end
        
        captitle=sprintf('%s',calprefix);
	    capfold=sprintf('%s.fat_band',calprefix);
	    folderc=sprintf('%s/',outdir);
	    
        if overwrite
            capnew=sprintf('%s/%s',outdir,capfold);
            [~, ~, ~] = mkdir(capnew);
        else
            if exist([folderc,capfold], 'dir')
                num=1;
                while exist([folderc,capfold, '_', num2str(num)], 'dir')
                    num=num+1;
                end
                new_fold=[capfold, '_', num2str(num)];
                capnew=sprintf('%s/%s',outdir,new_fold);
                [~, ~, ~] = mkdir(capnew);
            else
                capnew=sprintf('%s/%s',outdir,capfold);
                [~, ~, ~] = mkdir(capnew);
            end	
        end
        
        folderc=sprintf('%s/',capnew);
        
        if overwrite
            new_filename = [captitle, '.txt']; 
        else
            if exist([folderc,captitle,'.txt'], 'file')
                number = 1; 
                while exist([folderc,captitle, '_', num2str(number), '.txt'], 'file')
                    number = number + 1; 
                end
                new_filename = [captitle, '_', num2str(number), '.txt'];
            else
                new_filename = [captitle, '.txt']; 
            end
        end
	    
        path=sprintf('%s/%s',capnew,new_filename);
        outdir=capnew;
        band_source=sprintf("%s/%s",inputdir,new_filename_band_in);
        copyfile(band_source, outdir);
        switch nspin_type
            case "no_soc"
                [numB,xstt,ystt,wstt]=fat_band_no_soc(prefix,path,inputdir,outdir,Efermi);
            case "soc"
                [numB,xstt,ystt,wstt]=fat_band_soc(prefix,path,inputdir,outdir,Efermi);
            case "spin_up_down"
                [numB,xstt,ystt,wstt]=fat_band_soc_spin_up_down(prefix,path,inputdir,outdir,Efermi);
        end
    
    
        function[numBjy,xstjy,ystjy,wstjy]=fat_band_no_soc(prefixjy,pathjy,inputdirjy,outdirjy,Efermijy)
	        tgrand=tic;
            tabn=15;
            fileID=fopen(pathjy,'w');    
            t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            %\033[1mBold Text\033[0m
        
	        fprintf(fileID,"%*sProgram Started on : %s\n",tabn,'',t);
            fprintf(fileID,"\n%*s*****************Author Information*****************\n\n",tabn,'');
	        fprintf(fileID,"%*sMD. NILOY KHAN\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
	        fprintf("%*s<strong>Program Started on : %s</strong>\n",tabn,'',t);
            fprintf("\n%*s*****************Author Information*****************\n\n",tabn,'');
            fprintf("%*s[\bMD. NILOY KHAN]\b\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
            fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
            fprintf("%*s****************************************************\n\n",tabn,'');
            tstart1=tic;
            [States_matrix,atom_nam]=state_info_wo_soc(prefixjy,inputdirjy);
            % atoms=unique(atom_nam);
            % States_matrix;
            % num_atom=max(States_matrix(:,2));
            tend1=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            fprintf("%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            %States_matrix(:,[4,5])=States_matrix(:,[5,4]);
            tstart1=tic;
            if Efermijy
                Fermi=get_Fermi_energy(prefixjy,inputdirjy);
            else
                Fermi=0;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            
            tstart1=tic;
            num_state=max(States_matrix(:,1));
            [psi_no,psi_wgt,band,kpoints]=awk_proj(prefixjy,inputdirjy,num_state);
            [kp,nbnd]=size(band);
            band=band-Fermi;
            tend3=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            fprintf("%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            
            tstart1=tic;
            numBjy=cell(length(atom_nam),1);
            for e=1:length(atom_nam)
                %new_cap=sprintf('%.3d. %s_%.2d_%d_%d_%d_% 0.2f',States_matrix(e,1),atom_nam{e},States_matrix(e,2),States_matrix(e,3),States_matrix(e,4),States_matrix(e,5),States_matrix(e,6));              
                new_cap=sprintf('state #  %.3d: atom   %.2d (%s ), wfc  %d (l=%d m= %d)',States_matrix(e,1),States_matrix(e,2),atom_nam{e},States_matrix(e,3),States_matrix(e,4),States_matrix(e,5));
                numBjy{e,1}=new_cap;         
            end
            tend5=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            fprintf("%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            
             tstart1=tic;
             psii_noo=cell(kp,1);
             psii_wgtt=cell(kp,1);
             for pp=1:kp
                 Awe=psi_no{pp};
                 Awe(cellfun(@isempty, Awe)) = {0};
                 Wei=psi_wgt{pp};
                 Wei(cellfun(@isempty, Wei)) = {0};
                 psii_noo{pp}=Awe;
                 psii_wgtt{pp}=Wei;
             end
        
            psi_wgt=psii_wgtt;
            % psi_no=psii_noo;
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
            
            
            tstart1=tic;
            xstjy=zeros(kp,nbnd,length(numBjy));
            ystjy=zeros(kp,nbnd,length(numBjy));
            wstjy=zeros(kp,nbnd,length(numBjy));
        
            for on=1:length(numBjy)
                ystjy(:,:,on)=band;
                xpos=repmat(kpoints',1,nbnd);
                xstjy(:,:,on)=xpos;
            end
            for p=1:kp
                weight=cell2mat(psi_wgt{p});
                wstjy(p,:,:)=weight;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
        
        
            fprintf(fileID,"%*sWriting Output to \n%*s%s/\n\n",tabn,'',tabn,'',outdirjy);
            fprintf("%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',outdirjy);
            tstart1=tic;
            cap_numB=sprintf('%s_state_name.txt',prefixjy);
            caption_numB=sprintf('%s/%s',outdirjy,cap_numB);
            writecell(numBjy,caption_numB)
            cap_xst=sprintf('%s_xposition.MAT',prefixjy);
            caption_xst=sprintf('%s/%s',outdirjy,cap_xst);
            save(caption_xst,"xstjy","-v7.3","-nocompression");
            cap_yst=sprintf('%s_yposition.MAT',prefixjy);
            caption_yst=sprintf('%s/%s',outdirjy,cap_yst);
            save(caption_yst,"ystjy","-v7.3","-nocompression");
            cap_wst=sprintf('%s_weight.MAT',prefixjy);
            caption_wst=sprintf('%s/%s',outdirjy,cap_wst);
            save(caption_wst,"wstjy","-v7.3","-nocompression");
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
		        
            tgraend=toc(tgrand);
            fprintf(fileID,"%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s*************************THE END***************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n\n",tabn,'');
	        
            fprintf("%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s*************************THE END***************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n\n",tabn,'');
            t=datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            fprintf(fileID,"%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fprintf("%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fclose(fileID);
        end
    
        function[numBzy,xstzy,ystzy,wstzy]=fat_band_soc(prefixzy,pathzy,inputdirzy,outdirzy,Efermizy)
	        tgrand=tic;
        
            tabn=15;
            fileID=fopen(pathzy,'w');    
            t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            %\033[1mBold Text\033[0m
        
	        fprintf(fileID,"%*sProgram Started on : %s\n",tabn,'',t);
            fprintf(fileID,"\n%*s*****************Author Information*****************\n\n",tabn,'');
	        fprintf(fileID,"%*sMD. NILOY KHAN\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
	        fprintf("%*s<strong>Program Started on : %s</strong>\n",tabn,'',t);
            fprintf("\n%*s*****************Author Information*****************\n\n",tabn,'');
            fprintf("%*s[\bMD. NILOY KHAN]\b\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
            fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
            fprintf("%*s****************************************************\n\n",tabn,'');
            tstart1=tic;
            [States_matrix,atom_nam]=state_info_w_soc(prefixzy,inputdirzy);
            % atoms=unique(atom_nam);
            % num_atom=max(States_matrix(:,2));
            tend1=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            fprintf("%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            %States_matrix(:,[4,5])=States_matrix(:,[5,4]);
            tstart1=tic;
            if Efermizy
                Fermi=get_Fermi_energy(prefixzy,inputdirzy);
            else
                Fermi=0;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            
            tstart1=tic;
            num_state=max(States_matrix(:,1));
            [psi_no,psi_wgt,band,kpoints]=awk_proj(prefixzy,inputdirzy,num_state);
            [kp,nbnd]=size(band);
            band=band-Fermi;
            tend3=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            fprintf("%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            
            tstart1=tic;
            numBzy=cell(length(atom_nam),1);
            for e=1:length(atom_nam)
                %new_cap=sprintf('%.3d. %s_%.2d_%d_%d_%d_% 0.2f',States_matrix(e,1),atom_nam{e},States_matrix(e,2),States_matrix(e,3),States_matrix(e,4),States_matrix(e,5),States_matrix(e,6));              
                new_cap=sprintf('state #  %.3d: atom   %.2d (%s ), wfc  %d (j=%0.1f l= %d m_j=% 0.1f)',States_matrix(e,1),States_matrix(e,2),atom_nam{e},States_matrix(e,3),States_matrix(e,4),States_matrix(e,5),States_matrix(e,6));
                numBzy{e,1}=new_cap;         
            end
            tend5=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            fprintf("%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            
            tstart1=tic;
            psii_noo=cell(kp,1);
            psii_wgtt=cell(kp,1);
             for pp=1:kp
                 Awe=psi_no{pp};
                 Awe(cellfun(@isempty, Awe)) = {0};
                 Wei=psi_wgt{pp};
                 Wei(cellfun(@isempty, Wei)) = {0};
                 psii_noo{pp}=Awe;
                 psii_wgtt{pp}=Wei;
             end
        
            psi_wgt=psii_wgtt;
            % psi_no=psii_noo;
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
            
            tstart1=tic;
            xstzy=zeros(kp,nbnd,length(numBzy));
            ystzy=zeros(kp,nbnd,length(numBzy));
            wstzy=zeros(kp,nbnd,length(numBzy));
            
            for on=1:length(numBzy)
                ystzy(:,:,on)=band;
                xpos=repmat(kpoints',1,nbnd);
                xstzy(:,:,on)=xpos;
            end
            for p=1:kp
                weight=cell2mat(psi_wgt{p});
                wstzy(p,:,:)=weight;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
        
        
            fprintf(fileID,"%*sWriting Output to \n%*s%s/\n\n",tabn,'',tabn,'',outdirzy);
            fprintf("%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',outdirzy);
            tstart1=tic;
            cap_numB=sprintf('%s_state_name.txt',prefixzy);
            caption_numB=sprintf('%s/%s',outdirzy,cap_numB);
            writecell(numBzy,caption_numB)
            cap_xst=sprintf('%s_xposition.MAT',prefixzy);
            caption_xst=sprintf('%s/%s',outdirzy,cap_xst);
            save(caption_xst,"xstzy","-v7.3","-nocompression");
            cap_yst=sprintf('%s_yposition.MAT',prefixzy);
            caption_yst=sprintf('%s/%s',outdirzy,cap_yst);
            save(caption_yst,"ystzy","-v7.3","-nocompression");
            cap_wst=sprintf('%s_weight.MAT',prefixzy);
            caption_wst=sprintf('%s/%s',outdirzy,cap_wst);
            save(caption_wst,"wstzy","-v7.3","-nocompression");
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
		        
            tgraend=toc(tgrand);
            fprintf(fileID,"%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s*************************THE END***************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n\n",tabn,'');
	        
            fprintf("%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s*************************THE END***************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n\n",tabn,'');
            t=datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            fprintf(fileID,"%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fprintf("%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fclose(fileID);
        end
    
        function[numBvy,xstvy,ystvy,wstvy]=fat_band_soc_spin_up_down(prefixvy,pathvy,inputdirvy,outdirvy,Efermivy)
	        tgrand=tic;
            tabn=15;
	        
            fileID=fopen(pathvy,'w');    
            t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            %\033[1mBold Text\033[0m
        
	        fprintf(fileID,"%*sProgram Started on : %s\n",tabn,'',t);
            fprintf(fileID,"\n%*s*****************Author Information*****************\n\n",tabn,'');
	        fprintf(fileID,"%*sMD. NILOY KHAN\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
	        fprintf("%*s<strong>Program Started on : %s</strong>\n",tabn,'',t);
            fprintf("\n%*s*****************Author Information*****************\n\n",tabn,'');
            fprintf("%*s[\bMD. NILOY KHAN]\b\n%*sMSc in EEE, BUET(on going)\n%*sEmail: khanniloy534@gmail.com\n%*sPhone: +8801619141513\n",tabn,'',tabn,'',tabn,'',tabn,'');
            
            fprintf(fileID,"%*s****************************************************\n\n",tabn,'');
            fprintf("%*s****************************************************\n\n",tabn,'');
            tstart1=tic;
            [States_matrix,atom_nam]=state_info_w_soc_spin_up_down(prefixvy,inputdirvy);
            % atoms=unique(atom_nam);
            % num_atom=max(States_matrix(:,2));
            tend1=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            fprintf("%*sTime spend for extracting State Information : %f seconds\n",tabn,'',tend1);
            %States_matrix(:,[4,5])=States_matrix(:,[5,4]);
            tstart1=tic;
            if Efermivy
                Fermi=get_Fermi_energy(prefixvy,inputdirvy);
            else
                Fermi=0;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for extracting Band Data : %f seconds\n",tabn,'',tend2);
            
            tstart1=tic;
            num_state=max(States_matrix(:,1));
            [psi_no,psi_wgt,band,kpoints]=awk_proj(prefixvy,inputdirvy,num_state);
            [kp,nbnd]=size(band);
            band=band-Fermi;
            tend3=toc(tstart1);
            fprintf(fileID,"%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            fprintf("%*sTime spend for extracting Wave Function : %f seconds\n",tabn,'',tend3);
            
            tstart1=tic;
            numBvy=cell(length(atom_nam),1);
            for e=1:length(atom_nam)
                %new_cap=sprintf('%.3d. %s_%.2d_%d_%d_%d_% 0.2f',States_matrix(e,1),atom_nam{e},States_matrix(e,2),States_matrix(e,3),States_matrix(e,4),States_matrix(e,5),States_matrix(e,6));              
                new_cap=sprintf('state #  %.3d: atom   %.2d (%s ), wfc  %d (l=%d m= %d s_z=% 0.1f)',States_matrix(e,1),States_matrix(e,2),atom_nam{e},States_matrix(e,3),States_matrix(e,4),States_matrix(e,5),States_matrix(e,6));
                numBvy{e,1}=new_cap;         
            end
            tend5=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            fprintf("%*sTime spend for generating distinct variable : %f seconds\n",tabn,'',tend5);
            
            tstart1=tic;
            psii_noo=cell(kp,1);
            psii_wgtt=cell(kp,1);
            for pp=1:kp
                Awe=psi_no{pp};
                Awe(cellfun(@isempty, Awe)) = {0};
                Wei=psi_wgt{pp};
                Wei(cellfun(@isempty, Wei)) = {0};
                psii_noo{pp}=Awe;
                psii_wgtt{pp}=Wei;
            end
        
            psi_wgt=psii_wgtt;
            % psi_no=psii_noo;
            tend2=toc(tstart1);
           fprintf(fileID,"%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
           fprintf("%*sTime spend for replacing empty weight cell with 0 : %f seconds\n",tabn,'',tend2);
            
            tstart1=tic;
            xstvy=zeros(kp,nbnd,length(numBvy));
            ystvy=zeros(kp,nbnd,length(numBvy));
            wstvy=zeros(kp,nbnd,length(numBvy));
            
            for on=1:length(numBvy)
                ystvy(:,:,on)=band;
                xpos=repmat(kpoints',1,nbnd);
                xstvy(:,:,on)=xpos;
            end
            for p=1:kp
                weight=cell2mat(psi_wgt{p});
                wstvy(p,:,:)=weight;
            end
            
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for generating weight as a function of kpoints, band & states : %f seconds\n",tabn,'',tend2);
        
        
            fprintf(fileID,"%*sWriting Output to \n%*s%s/\n\n",tabn,'',tabn,'',outdirvy);
            fprintf("%*sWriting Output to \n%*s%s\n\n",tabn,'',tabn,'',outdirvy);
            tstart1=tic;
            cap_numB=sprintf('%s_state_name.txt',prefixvy);
            caption_numB=sprintf('%s/%s',outdirvy,cap_numB);
            writecell(numBvy,caption_numB)
            cap_xst=sprintf('%s_xposition.MAT',prefixvy);
            caption_xst=sprintf('%s/%s',outdirvy,cap_xst);
            save(caption_xst,"xstvy","-v7.3","-nocompression");
            cap_yst=sprintf('%s_yposition.MAT',prefixvy);
            caption_yst=sprintf('%s/%s',outdirvy,cap_yst);
            save(caption_yst,"ystvy","-v7.3","-nocompression");
            cap_wst=sprintf('%s_weight.MAT',prefixvy);
            caption_wst=sprintf('%s/%s',outdirvy,cap_wst);
            save(caption_wst,"wstvy","-v7.3","-nocompression");
            tend2=toc(tstart1);
            fprintf(fileID,"%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
            fprintf("%*sTime spend for writing weight,kpoint,band : %f seconds\n",tabn,'',tend2);
		        
            tgraend=toc(tgrand);
            fprintf(fileID,"%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s*************************THE END***************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n",tabn,'');
            fprintf(fileID,"%*s***********************************************************\n\n",tabn,'');
	        
            fprintf("%*sExecution Time of the Programs                : %.2f seconds\n\n",tabn,'',tgraend);
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s*************************THE END***************************\n",tabn,'');
            fprintf("%*s***********************************************************\n",tabn,'');
            fprintf("%*s***********************************************************\n\n",tabn,'');
            t=datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss');
            
            fprintf(fileID,"%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fprintf("%*sProgram Terminated on                 : %s\n",tabn,'',t);
            fclose(fileID);
        end
        
        function [psi,psi_weight,band_new,kpoints]=awk_proj(prefixny,inputdirny,num_stateny)
        % clc
        % clear
        % close all
        %     inputdir="./";
        %     prefix="Graphene_soc";
        %     num_state=16;
            new_filenameny=sprintf('%s.pdos.out',prefixny);
            capname44=sprintf('%s/%s',inputdirny,new_filenameny);
        
            fidny = fopen(capname44, 'r');
            data=textscan(fidny,'%s','whitespace','+','whitespace','*[#]+');
            data=data{1,1};
            fclose(fidny);
            % search for the k-points
            k_indexx = find(contains(data, 'k ='));
            k_index =[k_indexx;length(data)];
            % kpointss=data(k_indexx);
            kplen=length(k_indexx);
        
            eww_index = contains(data, '==== e(');
            band=data(eww_index);
            band_new = cellfun(@(x) sscanf(x, '==== e(%*d) = %f eV ===='), band);
            nnbnd=length(band_new)/kplen;
            % psi2_index = contains(data, '|psi|^2 =');
            % psi2=data(psi2_index);
            % psi2_new = cellfun(@(x) sscanf(x, '|psi|^2 = %d'), psi2);
            
            
            psi_weight = cell(kplen, 1);
            psi=cell(kplen,1);
            
            % Populate each cell of the outer array with an inner cell array
            for iny = 1:kplen
                for jny = 1:1
                    psi_weight{iny,jny} = cell(nnbnd, num_stateny);
                end
            end
            % psi_wwt=psi_weight;
            
            for k=1:kplen
                ew_index = find(contains(data(k_index(k):k_index(k+1)), '==== e('));
                %data(ew_index)
                psi_index = find(contains(data(k_index(k):k_index(k+1)), '|psi|^2 ='));
                d=0;
                for ew=1:nnbnd
                    p=0;
                    d=d+1;
                    for rg=k_index(k)+ew_index(ew):2:k_index(k)+psi_index(ew)-2
                        p=p+1;
			            % mystrr{d,p}=data{rg};
			            mystr=data{rg};
			            % mystrr2{d,p}=data{rg+1};
			            mystr2=data{rg+1};
                        value1=str2double(mystr2);
                        if isnan(value1)
                           psi{k,1}{ew,p}=[];
                        else
                           psi{k,1}{ew,p}=value1;
                        end
        				            
                        if  p==1
                            stind=strfind(mystr,'=');
                            value=str2double(mystr(stind+2:end));
                            if isnan(value)
	                            psi_weight{k,1}{ew,p}=0;
                            else
	                            psi_weight{k,1}{ew,value1}=value;
                            end 
                        else
                            value=str2double(mystr);
                            if isnan(value)
                                psi_weight{k,1}{ew,p}=0;
                            else
	                            psi_weight{k,1}{ew,value1}=value;
                            end
                        end
                    end
                end
            end
            kpoints=linspace(1,kplen,kplen);
            band_new=reshape(band_new,nnbnd,kplen);
            % psi2_new=reshape(psi2_new,nnbnd,kplen);
            band_new=band_new';
            % psi2_new=psi2_new';
        
            % figure
            % plot(band_new)
        end
        
        function[Fermi]=get_Fermi_energy(prefixey,inputdirey)
            %%able to generate bands data from only nscf output with low verbosity & band +scf output with high/low verbosity
            %%required scf output for fermi energy info if del ='band'
            %%%%% name of output file
            %prefix='Na3Bi_ml';
            %del='band';
        
            new_filename6 = sprintf('%s.scf.out',prefixey);
            caption6=sprintf('%s/%s',inputdirey,new_filename6);
            %%%%%%%%
            
            ffid=fopen(caption6);
            Caap=textscan(ffid,'%s','delimiter','\n');
            index5=find(contains(Caap{1},'the Fermi energy is'),1);
            index6=find(contains(Caap{1},'highest occupied, lowest unoccupied level (ev):'),1);
            
            if ~isempty(index5)
                frewind(ffid)
                Cap=textscan(ffid,'%s','headerlines',index5-1,'whitespace','the Fermi Energy is','collectoutput',true);
                fermi=Cap{1,1}{1,1};
                Fermi=str2double(fermi);
            elseif ~isempty(index6)
                frewind(ffid)
                Cap=textscan(ffid,'%f%f','headerlines',index6-1,'whitespace','highest occupied, lowest unoccupied level (ev):','collectoutput',true);
                fermi=cell2mat(Cap);
                Fermi=fermi(1)+(fermi(2)-fermi(1))/2;
            end
            
            frewind(ffid)  
            
            fclose(ffid);
        end
        
        function[state_matrix,atom_nam]=state_info_wo_soc(prefixty,inputdirty)
            % clc
            % clear 
            % close all
            % prefix='Bi4Br4_relax_final';
            % inputdir="./";
            new_filenamety=sprintf('%s.pdos.out',prefixty);
            capname44=sprintf('%s/%s',inputdirty,new_filenamety);
            % Open the text file for reading
            fidty = fopen(capname44, 'r');
        
            num_state=0;
            while ~feof(fidty)
                line = fgetl(fidty);
                if ~isempty(line)
                   match_for_ind=regexp(line,'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+)\s+m=\s*\-?(\d+)\)','tokens');
                    if ~isempty(match_for_ind)
                        num_state=num_state+1;
                    end
                end
            end
            frewind(fidty);
            
            % Initialize empty matrices to store the state information
            state_nums = zeros(1,num_state);
            atom_nums = zeros(1,num_state);
            atom_nam=cell(1,num_state);
            wfc_nums = zeros(1,num_state);
            l_vals = zeros(1,num_state);
            m_vals = zeros(1,num_state);
            %mj_vals = [];
            
            % Read each line of the text file and extract the state information
            k=0;
            while ~feof(fidty)
                line = fgetl(fidty);
                if ~isempty(line)
                    % Use regular expressions to extract the state information from the line
                    %match = regexp(line, '(\d+), (\d+), (\d+), (\d+), ([\d\.]+), ([\-\d\.]+);', 'tokens')
                    %match = regexp(line, 'state # (\d+): atom (\d+) \(.*\), wfc (\d+) \(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens')
                    %match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\(.*\),\s+wfc\s+(\d+)\s+\(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens');
                   % match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+) m=([\d\.]+)\)', 'tokens');
                   match=regexp(line,'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+)\s+m=\s*\-?(\d+)\)','tokens');
                    if ~isempty(match) %|| ~isempty(match1)
                        % str2double(match{1}{1})
                        k=k+1;
                        state_nums(k) = str2double(match{1}{1});
                        atom_nam{k}   =match{1}{3};
                        atom_nums(k) = str2double(match{1}{2});
                        wfc_nums(k) = str2double(match{1}{4});
                        l_vals(k) = str2double(match{1}{5});
                        m_vals(k) = str2double(match{1}{6});
                       % mj_vals(end+1) = str2double(match{1}{7});
                    end
                end
            end
            state_matrix = zeros(length(state_nums), 5);
            state_matrix(:,1) = state_nums';
            state_matrix(:,2) = atom_nums';
            state_matrix(:,3) = wfc_nums';
            state_matrix(:,4) = l_vals';
            state_matrix(:,5) = m_vals';
            %state_matrix(:,6) = mj_vals';
            %writematrix(state_matrix,'States_information_matrix.txt','delimiter','tab')
            
            % Close the file
            fclose(fidty);
        end
    
        function[state_matrixay,atom_namay]=state_info_w_soc(prefixay,inputdiray)
            %clc
            %clear all
            %close all
            %prefix='Bi4Br4_ml_soc'
            new_filenameay=sprintf('%s.pdos.out',prefixay);
            capname44=sprintf('%s/%s',inputdiray,new_filenameay);
            % Open the text file for reading
            fiday = fopen(capname44, 'r');
            
            num_state=0;
            while ~feof(fiday)
                line = fgetl(fiday);
                if ~isempty(line)
                   match_for_ind =  regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(j=(\d+\.\d+)\s+l=(\d+)\s+m_j=\s*([-]?\d+\.\d+)\)','tokens');
                    if ~isempty(match_for_ind)
                        num_state=num_state+1;
                    end
                end
            end
            frewind(fiday);
            
            % Initialize empty matrices to store the state information
            state_nums = zeros(1,num_state);
            atom_nums = zeros(1,num_state);
            atom_namay=cell(1,num_state);
            wfc_nums =zeros(1,num_state);
            l_vals = zeros(1,num_state);
            j_vals = zeros(1,num_state);
            mj_vals = zeros(1,num_state);
            
            k=0;
            % Read each line of the text file and extract the state information
            while ~feof(fiday)
                line = fgetl(fiday);
                if ~isempty(line)
                    % Use regular expressions to extract the state information from the line
                    %match = regexp(line, '(\d+), (\d+), (\d+), (\d+), ([\d\.]+), ([\-\d\.]+);', 'tokens')
                    %match = regexp(line, 'state # (\d+): atom (\d+) \(.*\), wfc (\d+) \(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens')
                    %match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\(.*\),\s+wfc\s+(\d+)\s+\(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens');
                    %match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+) j=([\d\.]+) m_j\s*=\s*([\-\d\.]+)\)', 'tokens');
                    match =  regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(j=(\d+\.\d+)\s+l=(\d+)\s+m_j=\s*([-]?\d+\.\d+)\)','tokens');
                    if ~isempty(match) %|| ~isempty(match1)
                        k=k+1;
                        state_nums(k) = str2double(match{1}{1});
                        atom_namay{k}   =match{1}{3};
                        atom_nums(k) = str2double(match{1}{2});
                        wfc_nums(k) = str2double(match{1}{4});
                        l_vals(k) = str2double(match{1}{5});
                        j_vals(k) = str2double(match{1}{6});
                        mj_vals(k) = str2double(match{1}{7});
                    end
                end
            end
            state_matrixay = zeros(length(state_nums), 6);
            state_matrixay(:,1) = state_nums';
            state_matrixay(:,2) = atom_nums';
            state_matrixay(:,3) = wfc_nums';
            state_matrixay(:,4) = l_vals';
            state_matrixay(:,5) = j_vals';
            state_matrixay(:,6) = mj_vals';
            %writematrix(state_matrix,'States_information_matrix.txt','delimiter','tab')
            
            % Close the file
            fclose(fiday);
        end
    
        function[state_matrixmy,atom_nammy]=state_info_w_soc_spin_up_down(prefixmy,inputdirmy)
            %clc
            %clear all
            %close all
            %prefix='Bi4Br4_ml_soc'
            new_filenamemy=sprintf('%s.pdos.out',prefixmy);
            capname44=sprintf('%s/%s',inputdirmy,new_filenamemy);
            % Open the text file for reading
            fidmy = fopen(capname44, 'r');
            num_state=0;
            while ~feof(fidmy)
                line = fgetl(fidmy);
                if ~isempty(line)
                   match_for_ind =  regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+)\s+m=( \d+)\s+s_z=\s*([-]?\d+\.\d+)\)','tokens');
                    if ~isempty(match_for_ind)
                        num_state=num_state+1;
                    end
                end
            end
            frewind(fidmy);
            % Initialize empty matrices to store the state information
            state_nums = zeros(1,num_state);
            atom_nums = zeros(1,num_state);
            atom_nammy=cell(1,num_state);
            wfc_nums = zeros(1,num_state);
            l_vals = zeros(1,num_state);
            j_vals = zeros(1,num_state);
            mj_vals = zeros(1,num_state);
            
            % Read each line of the text file and extract the state information
            k=0;
            while ~feof(fidmy)
                line = fgetl(fidmy);
                if ~isempty(line)
                    % Use regular expressions to extract the state information from the line
                    %match = regexp(line, '(\d+), (\d+), (\d+), (\d+), ([\d\.]+), ([\-\d\.]+);', 'tokens')
                    %match = regexp(line, 'state # (\d+): atom (\d+) \(.*\), wfc (\d+) \(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens')
                    %match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\(.*\),\s+wfc\s+(\d+)\s+\(l=(\d+) j=([\d\.]+) m_j=([\-\d\.]+)\)', 'tokens');
                    %match = regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+) j=([\d\.]+) m_j\s*=\s*([\-\d\.]+)\)', 'tokens');
                    match =  regexp(line, 'state #\s+(\d+):\s+atom\s+(\d+)\s+\((\w+)\s*\),\s+wfc\s+(\d+)\s+\(l=(\d+)\s+m=( \d+)\s+s_z=\s*([-]?\d+\.\d+)\)','tokens');
                    if ~isempty(match) %|| ~isempty(match1)
                        k=k+1;
                        state_nums(k) = str2double(match{1}{1});
                        atom_nammy{k}   =match{1}{3};
                        atom_nums(k) = str2double(match{1}{2});
                        wfc_nums(k) = str2double(match{1}{4});
                        l_vals(k) = str2double(match{1}{5});
                        j_vals(k) = str2double(match{1}{6});
                        mj_vals(k) = str2double(match{1}{7});
                    end
                end
            end
            state_matrixmy = zeros(length(state_nums), 6);
            state_matrixmy(:,1) = state_nums';
            state_matrixmy(:,2) = atom_nums';
            state_matrixmy(:,3) = wfc_nums';
            state_matrixmy(:,4) = l_vals';
            state_matrixmy(:,5) = j_vals';
            state_matrixmy(:,6) = mj_vals';
            %writematrix(state_matrix,'States_information_matrix.txt','delimiter','tab')
            
            % Close the file
            fclose(fidmy);
        end
    end
    
    
    function [ppp]=fat_band_plot(varargin)  
        newvararginnew = varargin;
        varargin = cell(1, 16); % Initialize with 34 empty cells
    
        for j = 1:length(newvararginnew)
            if isempty(newvararginnew{j})
                % If input is an empty array, store empty array
                varargin{j} = [];
            else
                % If input is zero or any other non-empty value, store it as is
                varargin{j} = newvararginnew{j};
            end
        end
        
        if nargin<17
            if isempty(varargin{1})
                error('prefix is required');
            else
                prefix=varargin{1};
                cap_xst=sprintf('%s_xposition.MAT',prefix);
            end
            if isempty(varargin{2})
                error('Input directory is required');
            else
                inputdir=varargin{2};
                caption_xst=sprintf('%s/%s',inputdir,cap_xst);
                xxxst=cell2mat(struct2cell(load(caption_xst)));
            end
            if isempty(varargin{3})
                error("must define want states in matrix form.\n" + ...
                    " each row represent a specific type of states.\n" + ...
                    "zero padding needed if want_states in each row is not equal.");
            else
                Want_state=varargin{3};
                [Wantr,Wantc]=size(Want_state);
            end
            if isempty(varargin{4})
                outdir='.';
            else
                outdir=varargin{4};
            end
    
            if isempty(varargin{5})
                band_limit=[1 length(xxxst(1,:,1))];
                bdstart=band_limit(1);
                bdend=band_limit(2);
                xst=xxxst(:,bdstart:bdend,:);
                cap_yst=sprintf('%s_yposition.MAT',prefix);
                caption_yst=sprintf('%s/%s',inputdir,cap_yst);
                yyyst=cell2mat(struct2cell(load(caption_yst)));
                yst=yyyst(:,bdstart:bdend,:);
                
                cap_wst=sprintf('%s_weight.MAT',prefix);
                caption_wst=sprintf('%s/%s',inputdir,cap_wst);
                wwwst=cell2mat(struct2cell(load(caption_wst)));
                wst=wwwst(:,bdstart:bdend,:);
            else
                band_limit=varargin{5};
                if band_limit(1)==0
                    bdstart=1;
                    bdend=band_limit(2);
                elseif band_limit(2)==0
                    bdstart=band_limit(1);
                    bdend=length(xxxst(1,:,1));
                elseif band_limit(1)==0 && band_limit(2)==0
                    bdstart=1;
                    bdend=length(xxxst(1,:,1));
                else
                    bdstart=band_limit(1);
                    bdend=band_limit(2);
                end
                xst=xxxst(:,bdstart:bdend,:);
                cap_yst=sprintf('%s_yposition.MAT',prefix);
                caption_yst=sprintf('%s/%s',inputdir,cap_yst);
                yyyst=cell2mat(struct2cell(load(caption_yst)));
                yst=yyyst(:,bdstart:bdend,:);
                
                cap_wst=sprintf('%s_weight.MAT',prefix);
                caption_wst=sprintf('%s/%s',inputdir,cap_wst);
                wwwst=cell2mat(struct2cell(load(caption_wst)));
                wst=wwwst(:,bdstart:bdend,:);
            end
    
            if isempty(varargin{6})
                % hsv_mat=hsv(Wantr);
                hsv_mat=rand(Wantr,3);
            else
                hsv_mat=varargin{6};
            end
    
            if isempty(varargin{7})
                as_des_type='ascend';
            else
                as_des_type=varargin{7};
            end
            if isempty(varargin{8})
                xmin=min(xst,[], 'all');
                xmax=max(xst,[], 'all');
                ymin=min(yst,[], 'all');
                ymax=max(yst,[], 'all');
            else
                xylimit=varargin{8};
                if ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                end
            end
            if isempty(varargin{9})
                interpolation=0;
            else
                interpolation=1;
                interpoint=varargin{9};
            end
            if isempty(varargin{10})
                weight_multiplier=100;
            else
                weight_multiplier=varargin{10};
            end
            
    
            if isempty(varargin{11})
                justify_zero=0.01;
            else
                justify_zero=varargin{11};
            end
    
            if isempty(varargin{12})
                figure_spec={'Times','bold',12};
            else
                figure_spec=varargin{12};
            end
    
            if isempty(varargin{13})
                draw_fermi=1;
            else
                draw_fermi=varargin{13};
            end
    
            % if isempty(varargin{21})
            % 
            %     draw_arrow=varargin{21};
            %     xarrow=varargin{22};
            %     yarrow=varargin{23};
            % end
            
            if isempty(varargin{14})
                cap_legend=generateAlphabetCellArrayy(Wantr);
                cap_spec={'Location','bestoutside','Times','bold',10};
            else
                legend_details=varargin{14};
                if length(legend_details)==Wantr+5
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:5
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                elseif length(legend_details)==Wantr+4
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:4
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+3
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:3
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,4}='bold';
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+2
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:2
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,3}="Times";
                    cap_spec{1,4}='bold';
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+1 || length(legend_details)==Wantr
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec={'Location','bestoutside','Times','bold',10};
                end
            end
    
            if isempty(varargin{15}) 
                pbmod_rat=[1 1 1];
            else
	            pbmod_rat=varargin{15};
            end
            if isempty(varargin{16})   
                collect_fig="true";
                cap_fig=prefix;
            else
                var16=varargin{16};
                collect_fig=var16{1};
                switch collect_fig
                    case "true"
                        cap_fig=var16{2};
                    case "fasle"
                end
            end
        end
        
        [kpath,ktick,~,~,~]=ktick_extractorr(capname_band_in);
    
        figure 
        plot(xst(:,:,1),yst(:,:,1),'Color',[0.7 0.7 0.7])
    
        for d=1:length(xst(1,:,1))
            hold on
            if interpolation
                xstunsort=zeros(interpoint,Wantr);%[];
                ystunsort=zeros(interpoint,Wantr);%[];
                wstunsort=zeros(interpoint,Wantr);%[];
            else
                xstunsort=zeros(length(xst(:,1,1)),Wantr);%[];
                ystunsort=zeros(length(xst(:,1,1)),Wantr);%[];
                wstunsort=zeros(length(xst(:,1,1)),Wantr);%[];
            end

            
            for on=1:Wantr
                xxst=xst(:,d,on);
                yyst=yst(:,d,on);
                
                wwst=0;
                for lk=1:Wantc     
                    if Want_state(on,lk)~=0
                        wwst=wwst+wst(:,d,Want_state(on,lk));
                    end
                end
                wwst= wwst./length(nonzeros(Want_state(on,:)));
                
                if interpolation
                    xi = linspace(min(xxst), max(xxst), interpoint);
                    yi = interp1(xxst, yyst, xi);
                    wi = interp1(xxst, wwst, xi);
                    xstunsort(:,on)=xi';
                    ystunsort(:,on)=yi';
                    wstunsort(:,on)=wi';
                else
                    xstunsort(:,on)=xxst;
                    ystunsort(:,on)=yyst;
                    wstunsort(:,on)=wwst;
                end     
            end
        
            [wstsort, wstsortind] = sort(wstunsort, 2,as_des_type);
            pp=gobjects(Wantr,1);
            for col=1:Wantr
                col_ind=wstsortind(:,col);
                hsv_map=hsv_mat(col_ind,:);
                xpos=xstunsort(:,col);
                ypos=ystunsort(:,col);
                wpos=wstsort(:,col);
                pp(col)=scatter(xpos, ypos,weight_multiplier.*wpos+justify_zero,hsv_map,'filled','MarkerFaceAlpha', 0.9,'MarkerEdgeAlpha',0.9);
                %pp.Marker='_';
                %pp(col)=patch([xpos' nan],[ypos' nan],[wpos' nan],[wpos' nan], 'edgecolor', 'interp','LineWidth',5,'LineJoin','round');
            end
         
            xlim([xmin xmax])
            ylim([ymin ymax])
            
            if draw_fermi
                line([xmin xmax], [0, 0], 'Color', [0,0,0],'LineStyle','--','LineWidth',1);
                hold on
            end
            
            % if exist('draw_arrow','var') == 1
            %     [incr,~]=size(xarrow);
            %     for p=1:incr
            %         annotation('textarrow',xarrow(p,:),yarrow(p,:))
            %     end
            % end

            set(gca,'xtick',ktick,'xticklabel',kpath,'fontname',figure_spec{1},'fontweight',figure_spec{2},'fontsize',figure_spec{3});
    
            
            % if exist('xaxlabel','var') == 1
            %     xlabel(xlabel_cap{1},'FontName',figure_spec{1},'fontweight',figure_spec{2},'fontsize',str2double(figure_spec{3}))
            % end
            ylabel('Energy, E-E_F (eV)','FontName',figure_spec{1},'fontweight',figure_spec{2},'fontsize',figure_spec{3})
    
        end
    
        hold off
        box on
        ppp=pp;
        for h=1:Wantr
            %ppp(h).Marker='none';
            % ppp(h).YData=[1000];
            % ppp(h).XData=[1000];
            ppp(h).CData=hsv_mat(h,:);
            % ppp(h).CDataMode='manual';
         end
        legend(ppp,cap_legend,cap_spec{1},cap_spec{2},'FontName',cap_spec{3},'fontweight',cap_spec{4},'fontsize',cap_spec{5})
        legend('boxoff');
    
        pbaspect(pbmod_rat)
    
        switch collect_fig
            case 'true'
                cap_graphics = sprintf('%s/%s_fat_band.png',outdir,cap_fig);
                exportgraphics(gca,cap_graphics,'Resolution',1000) ;
            % exportgraphics(gca,'2d.eps','BackgroundColor','none','ContentType','vector','Resolution',3000)
        end
    end

    function[ppp]=fat_band_plot_col_var(varargin)  
        
        newvararginnew = varargin;
        varargin = cell(1, 16); % Initialize with 34 empty cells
    
        for j = 1:length(newvararginnew)
            if isempty(newvararginnew{j})
                % If input is an empty array, store empty array
                varargin{j} = [];
            else
                % If input is zero or any other non-empty value, store it as is
                varargin{j} = newvararginnew{j};
            end
        end
        
        if nargin<17
            if isempty(varargin{1})
                error('prefix is required');
            else
                prefix=varargin{1};
                cap_xst=sprintf('%s_xposition.MAT',prefix);
            end
            if isempty(varargin{2})
                error('Input directory is required');
            else
                inputdir=varargin{2};
                caption_xst=sprintf('%s/%s',inputdir,cap_xst);
                xxxst=cell2mat(struct2cell(load(caption_xst)));
            end
            if isempty(varargin{3})
                error("must define want states in matrix form.\n" + ...
                    " each row represent a specific type of states.\n" + ...
                    "zero padding needed if want_states in each row is not equal.");
            else
                Want_state=varargin{3};
                [Wantr,Wantc]=size(Want_state);
            end
            if isempty(varargin{4})
                outdir='.';
            else
                outdir=varargin{4};
            end
    
            if isempty(varargin{5})
                band_limit=[1 length(xxxst(1,:,1))];
                bdstart=band_limit(1);
                bdend=band_limit(2);
                xst=xxxst(:,bdstart:bdend,:);
                cap_yst=sprintf('%s_yposition.MAT',prefix);
                caption_yst=sprintf('%s/%s',inputdir,cap_yst);
                yyyst=cell2mat(struct2cell(load(caption_yst)));
                yst=yyyst(:,bdstart:bdend,:);
                
                cap_wst=sprintf('%s_weight.MAT',prefix);
                caption_wst=sprintf('%s/%s',inputdir,cap_wst);
                wwwst=cell2mat(struct2cell(load(caption_wst)));
                wst=wwwst(:,bdstart:bdend,:);
            else
                band_limit=varargin{5};
                if band_limit(1)==0
                    bdstart=1;
                    bdend=band_limit(2);
                elseif band_limit(2)==0
                    bdstart=band_limit(1);
                    bdend=length(xxxst(1,:,1));
                elseif band_limit(1)==0 && band_limit(2)==0
                    bdstart=1;
                    bdend=length(xxxst(1,:,1));
                else
                    bdstart=band_limit(1);
                    bdend=band_limit(2);
                end
                xst=xxxst(:,bdstart:bdend,:);
                cap_yst=sprintf('%s_yposition.MAT',prefix);
                caption_yst=sprintf('%s/%s',inputdir,cap_yst);
                yyyst=cell2mat(struct2cell(load(caption_yst)));
                yst=yyyst(:,bdstart:bdend,:);
                
                cap_wst=sprintf('%s_weight.MAT',prefix);
                caption_wst=sprintf('%s/%s',inputdir,cap_wst);
                wwwst=cell2mat(struct2cell(load(caption_wst)));
                wst=wwwst(:,bdstart:bdend,:);
            end
    
            if isempty(varargin{6})
                % hsv_mat=hsv(Wantr);
                hsv_mat=rand(Wantr,3);
            else
                hsv_mat=varargin{6};
            end
    
            if isempty(varargin{7})
                as_des_type='ascend';
            else
                as_des_type=varargin{7};
            end
            if isempty(varargin{8})
                xmin=min(xst,[], 'all');
                xmax=max(xst,[], 'all');
                ymin=min(yst,[], 'all');
                ymax=max(yst,[], 'all');
            else
                xylimit=varargin{8};
                if ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif ~isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=xylimit(1,1);%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && ~isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=xylimit(1,2);%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && ~isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=xylimit(2,1);
                    ymax=max(yst,[], 'all');
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && ~isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=xylimit(2,2);
                elseif isnan(xylimit(1,1)) && isnan(xylimit(1,2)) && isnan(xylimit(2,1)) && isnan(xylimit(2,2))
                    xmin=min(xst,[], 'all');%*min(xxxst,[], 'all')/length(xxxst(:,1,1));
                    xmax=max(xst,[], 'all');%*max(xxxst,[], 'all')/length(xxxst(:,1,1));
                    ymin=min(yst,[], 'all');
                    ymax=max(yst,[], 'all');
                end
            end
            if isempty(varargin{9})
                interpolation=0;
            else
                interpolation=1;
                interpoint=varargin{9};
            end
            if isempty(varargin{10})
                weight_multiplier=100;
            else
                weight_multiplier=varargin{10};
            end
            
    
            if isempty(varargin{11})
                justify_zero=0.01;
            else
                justify_zero=varargin{11};
            end
    
            if isempty(varargin{12})
                figure_spec={'Times','bold',12};
            else
                figure_spec=varargin{12};
            end
    
            if isempty(varargin{13})
                draw_fermi=1;
            else
                draw_fermi=varargin{13};
            end
    
            % if isempty(varargin{21})
            % 
            %     draw_arrow=varargin{21};
            %     xarrow=varargin{22};
            %     yarrow=varargin{23};
            % end
            
            if isempty(varargin{14})
                cap_legend=generateAlphabetCellArrayy(Wantr);
                cap_spec={'Location','bestoutside','Times','bold',10};
            else
                legend_details=varargin{14};
                if length(legend_details)==Wantr+5
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:5
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                elseif length(legend_details)==Wantr+4
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:4
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+3
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:3
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,4}='bold';
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+2
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec=cell(1,5);
                    for gj=1:2
                        cap_spec{1,gj}=legend_details{Wantr+gj};
                    end
                    cap_spec{1,3}="Times";
                    cap_spec{1,4}='bold';
                    cap_spec{1,5}=10;
                elseif length(legend_details)==Wantr+1 || length(legend_details)==Wantr
                    cap_legend=cell(1,Wantr);
                    for gj=1:Wantr
                        cap_legend{1,gj}=legend_details{gj};
                    end
                    cap_spec={'Location','bestoutside','Times','bold',10};
                end
            end
    
            if isempty(varargin{15}) 
                pbmod_rat=[1 1 1];
            else
                pbmod_rat=varargin{15};
            end
            if isempty(varargin{16})   
                collect_fig="true";
                cap_fig=prefix;
            else
                var16=varargin{16};
                collect_fig=var16{1};
                switch collect_fig
                    case "true"
                        cap_fig=var16{2};
                    case "fasle"
                end
            end
        end
        
        [kpath,ktick,~,~,~]=ktick_extractorr(capname_band_in);
        figure 
        plot(xst(:,:,1),yst(:,:,1),'Color',[0.7 0.7 0.7])
        
        for d=1:length(xst(1,:,1))
            hold on
            if interpolation
                xstunsort=zeros(interpoint,Wantr);%[];
                ystunsort=zeros(interpoint,Wantr);%[];
                wstunsort=zeros(interpoint,Wantr);%[];
            else
                xstunsort=zeros(length(xst(:,1,1)),Wantr);%[];
                ystunsort=zeros(length(xst(:,1,1)),Wantr);%[];
                wstunsort=zeros(length(xst(:,1,1)),Wantr);%[];
            end
            
            for on=1:Wantr
                xxst=xst(:,d,on);
                yyst=yst(:,d,on);
                wwst=0;
            
                for lk=1:Wantc     
                    if Want_state(on,lk)~=0
                        wwst=wwst+wst(:,d,Want_state(on,lk));
                    end
                end
                wwst= wwst./length(nonzeros(Want_state(on,:)));
            
                if interpolation
                    xi = linspace(min(xxst), max(xxst), interpoint);
                    yi = interp1(xxst, yyst, xi);
                    wi = interp1(xxst, wwst, xi);
                    xstunsort(:,on)=xi';
                    ystunsort(:,on)=yi';
                    wstunsort(:,on)=wi';
                else
                    xstunsort(:,on)=xxst;
                    ystunsort(:,on)=yyst;
                    wstunsort(:,on)=wwst;
                end     
            end
	        hsv_map_sum=0;
            [wstsort, wstsortind] = sort(wstunsort, 2,as_des_type);
	        wst_sum=sum(wstsort,2);
	        
	        for col=1:Wantr
		        col_ind=wstsortind(:,col);
                hsv_map=hsv_mat(col_ind,:);
                
		        wwgst=wstsort(:,col);
                hsv_store=zeros(length(wwgst),3);
		        for jk=1:length(wwgst)
			        wst_per=wwgst(jk)*1/wst_sum(jk);
			        hsv_store(jk,:)=wst_per*hsv_map(jk,:);
		        end
		        hsv_map_sum=hsv_map_sum+hsv_store;
	        end
            pp=gobjects(Wantr,1);
            for col=1:Wantr
                %col_ind=wstsortind(:,col);
                %hsv_map=hsv_mat(col_ind,:);
                xpos=xstunsort(:,col);
                ypos=ystunsort(:,col);
		        
                wpos=wstsort(:,col);
                pp(col)=scatter(xpos, ypos,weight_multiplier.*wpos+justify_zero,hsv_map_sum,'filled','MarkerFaceAlpha', 0.9,'MarkerEdgeAlpha',0.9);
                %pp.Marker='_';
                %pp(col)=patch([xpos' nan],[ypos' nan],[wpos' nan],[wpos' nan], 'edgecolor', 'interp','LineWidth',5,'LineJoin','round');
            end
        
        
            xlim([xmin xmax])
            ylim([ymin ymax])
            
            if draw_fermi
                line([xmin xmax], [0, 0], 'Color', [0,0,0],'LineStyle','--','LineWidth',1);
                hold on
            end
            
            % if exist('draw_arrow','var') == 1
            %     [incr,~]=size(xarrow);
            %     for p=1:incr
            %         annotation('textarrow',xarrow(p,:),yarrow(p,:))
            %     end
            % end
    
            set(gca,'xtick',ktick,'xticklabel',kpath,'fontname',figure_spec{1},'fontweight',figure_spec{2},'fontsize',figure_spec{3});
    
            
            % if exist('xaxlabel','var') == 1
            %     xlabel(xlabel_cap{1},'FontName',figure_spec{1},'fontweight',figure_spec{2},'fontsize',str2double(figure_spec{3}))
            % end
            ylabel('Energy, E-E_F (eV)','FontName',figure_spec{1},'fontweight',figure_spec{2},'fontsize',figure_spec{3})
    
        end
    
    % hold off
    % [hsvr,~]=size(hsv_mat);
    % for h=1:hsvr
    %     scaterp=pp(h);
    %     updatecdata=repmat(hsv_mat(h,:),length(ypos),1);
    %     scaterp.CData=updatecdata;
    %     pp(d)=scaterp;
    % end
    
        hold off
        box on
        ppp=pp;
        for h=1:Wantr
            %ppp(h).Marker='none';
            % ppp(h).YData=[1000];
            % ppp(h).XData=[1000];
            ppp(h).CData=hsv_mat(h,:);
            % ppp(h).CDataMode='manual';
         end
        legend(ppp,cap_legend,cap_spec{1},cap_spec{2},'FontName',cap_spec{3},'fontweight',cap_spec{4},'fontsize',cap_spec{5})
        legend('boxoff');
    
        pbaspect(pbmod_rat)

        if Wantr>1
            % Number of steps to interpolate between colors
            totalSteps = 256; % Number of interpolation steps for smooth gradient
            
            % Create the interpolated colormap
            interpColors = interp1(linspace(0, 1, size(hsv_mat, 1)), hsv_mat, linspace(0, 1, totalSteps), 'linear');
            
            % Apply the interpolated colormap
            colormap(interpColors);
            
            % Add color bar
            pppcol=colorbar;
            pppcol.Ticks=[];
        end
    
        switch collect_fig
            case 'true'
                cap_graphics = sprintf('%s/%s_col_var_fat_band.png',outdir,cap_fig);
                exportgraphics(gca,cap_graphics,'Resolution',1000) ;
            % exportgraphics(gca,'2d.eps','BackgroundColor','none','ContentType','vector','Resolution',3000)
        end
    end


    
    function varargout = ktick_extractorr(capname_band_in)
        fid = fopen(capname_band_in, 'r');
        data=textscan(fid,'%s','whitespace','+','whitespace','*[#]+');
        data=data{1,1};
        fclose(fid);
        nbnd_ind = find(contains(data, 'K_POINTS'));
        kpstring=strsplit(data{nbnd_ind});
        kpath_text = kpstring{2};
        kpath_text = strrep(kpath_text, '{', '');
        kpath_text = strrep(kpath_text, '}', '');
        kpath_text = strrep(kpath_text, '(', '');
        kpath_text = strrep(kpath_text, ')', '');
        kpath_text=lower(kpath_text);
    
        switch kpath_text
            case 'tpiba_b'
                Knumber = str2double(data{nbnd_ind+1});
            case 'crystal_b'
                Knumber = str2double(data{nbnd_ind+1});
            case 'crystal'
                Knumber = str2double(data{nbnd_ind+1});
            case 'automatic'
                mp_grid=str2double(strsplit(data{nbnd_ind+1}));
        end
    
    
    
        switch kpath_text
            case 'tpiba_b'
                kpathh=cell(1,Knumber);
                kp_pointss=zeros(1,Knumber);
            case 'crystal_b'
                kcord=zeros(Knumber,3);
                kp_pointss=zeros(1,Knumber);
            case 'crystal'
                kcord=zeros(Knumber,3);
            case 'automatic'
                kcord=zeros(mp_grid(1)*mp_grid(2)*mp_grid(3),3);
                % kp_pointss=zeros(1,mp_grid(1)*mp_grid(2)*mp_grid(3));
        end
        
        
        switch kpath_text
            case 'tpiba_b'
                jj=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    jj=jj+1;
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    switch components{1}
                        case 'gG'
                            kpathh{jj} = 'Γ';
                            kp_pointss(jj) = str2double(components{2});
                        otherwise
                            kpathh{jj} = components{1};
                            kp_pointss(jj) = str2double(components{2});
                    end
                end
                        
                kptick_rep= kp_pointss;
                kptick_rep(kptick_rep == 0) = 1;
                kpathtickk=zeros(1,Knumber);
                kpathtickk(1)=1;
                
                for g=1:Knumber-1
                    kpathtickk(g+1)=kpathtickk(g)+kptick_rep(g);
                end
            
                % Initialize the new array
                kpathtickk_array = kpathtickk;
                kpathtickk_org = kpathtickk;
            
                % Iterate through the array to adjust values
                % for i = 2:length(kpathtickk)
                %     diff = kpathtickk(i) - kpathtickk(i-1);
                %     if diff == 1 || diff == 2 || diff == 3
                %         % Ensure minimum separation of 2
                %         kpathtickk_array(i-1) = kpathtickk(i-1)-2;
                %         kpathtickk_array(i) = kpathtickk(i)+2;
                %         kpathtickk_org(i-1) = kpathtickk(i)- 2;
                %         kpathtickk_org(i) = kpathtickk(i)+ 2;
                % 
                %     else
                %         kpathtickk_array(i) = kpathtickk(i);
                %         kpathtickk_org(i) = kpathtickk(i);
                %     end
                % end
            case 'crystal_b'
                jj=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    jj=jj+1;
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    kcord(jj,1) = str2double(components{1});
                    kcord(jj,2) = str2double(components{2});
                    kcord(jj,3) = str2double(components{3});
                    kp_pointss(jj) = str2double(components{4});
                    if length(components)>4
                        temp_string=strrep(components{end},'!','');
                        if ~isempty(temp_string)
                            kpathh{jj}=temp_string;
                        end
                    end
    
                end
    
                kptick_rep= kp_pointss;
                kptick_rep(kptick_rep == 0) = 1;
                kpathtickk=zeros(1,Knumber);
                kpathtickk(1)=1;
                
                for g=1:Knumber-1
                    kpathtickk(g+1)=kpathtickk(g)+kptick_rep(g);
                end
            
                % Initialize the new array
                kpathtickk_array = kpathtickk;
                kpathtickk_org = kpathtickk;
            
                % Iterate through the array to adjust values
                % for i = 2:length(kpathtickk)
                %     diff = kpathtickk(i) - kpathtickk(i-1);
                %     if diff == 1 || diff == 2 || diff == 3
                %         % Ensure minimum separation of 2
                %         kpathtickk_array(i-1) = kpathtickk(i-1)-2;
                %         kpathtickk_array(i) = kpathtickk(i)+2;
                %         kpathtickk_org(i-1) = kpathtickk(i)- 2;
                %         kpathtickk_org(i) = kpathtickk(i)+ 2;
                % 
                %     else
                %         kpathtickk_array(i) = kpathtickk(i);
                %         kpathtickk_org(i) = kpathtickk(i);
                %     end
                % end
            case 'crystal'
                len_ktick=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    if length(components)>4
                        temp_string=strrep(components{end},'!','');
                        if ~isempty(temp_string)
                            len_ktick=len_ktick+1;
                        end
                    end
                end
                Ktick_index=zeros(1,len_ktick);
                jj=0;
                jjj=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    jj=jj+1;
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    kcord(jj,1) = str2double(components{1});
                    kcord(jj,2) = str2double(components{2});
                    kcord(jj,3) = str2double(components{3});
                    if length(components)>4
                        temp_string=strrep(components{end},'!','');
                        if ~isempty(temp_string)
                            jjj=jjj+1;
                            kpathh{jjj}=temp_string;
                            Ktick_index(jjj)=jj;
                        end
                    end
                end
                [kr,~]=size(kcord);
                kcord=2*pi.*kcord;
                if exist('Ktick_index','var') == 0
                    kpathtickk=[1 kr];   
                else
                    kpathtickk=Ktick_index;
                end
                kpathtickk_org=kpathtickk;
                kpathtickk_array=kpathtickk;
                kp_pointss=kr-1;
            case 'tpiba'
                len_ktick=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    if length(components)>4
                        temp_string=strrep(components{end},'!','');
                        if ~isempty(temp_string)
                            len_ktick=len_ktick+1;
                        end
                    end
                end
                Ktick_index=zeros(1,len_ktick);
                jj=0;
                jjj=0;
                for jvdg=nbnd_ind+2:nbnd_ind+Knumber+1
                    jj=jj+1;
                    linee=data{jvdg,1};
                    components = strsplit(linee);
                    kcord(jj,1) = str2double(components{1});
                    kcord(jj,2) = str2double(components{2});
                    kcord(jj,3) = str2double(components{3});
                    if length(components)>4
                        temp_string=strrep(components{end},'!','');
                        if ~isempty(temp_string)
                            jjj=jjj+1;
                            kpathh{jjj}=temp_string;
                            Ktick_index(jjj)=jj;
                        end
                    end
                end
                [kr,~]=size(kcord);
                if exist('Ktick_index','var') == 0
                    kpathtickk=[1 kr];   
                else
                    kpathtickk=Ktick_index;
                end
                kpathtickk_org=kpathtickk;
                kpathtickk_array=kpathtickk;
                kp_pointss=kr-1;
            case 'gamma'
                kcord=[0 0 0];
                kpathtickk=1;
                kpathtickk_org=kpathtickk;
                kpathtickk_array=kpathtickk;
                kp_pointss=1;
            case 'automatic'
                kx=mp_grid(1);
                ky=mp_grid(2);
                kz=mp_grid(3);
                kxx=linspace(0,1-1/kx,kx);
                kyy=linspace(0,1-1/ky,ky);
                kzz=linspace(0,1-1/kz,kz);
                [KX,KY,KZ]=meshgrid(kxx',kyy',kzz');
                KX=KX(:);
                KY=KY(:);
                KZ=KZ(:);
                kcrystal=[KX,KY,KZ];
                [kr,~]=size(kcrystal);
                kcord=2*pi.*kcrystal;
                kpathtickk=[1 kr];
                kpathtickk_org=kpathtickk;
                kpathtickk_array=kpathtickk;
                kp_pointss=kr-1;
        end
    
    
    
        switch kpath_text
            case 'tpiba_b'
                varargout=cell(1,5);
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
            case 'crystal_b'
                if exist('kpathh','var') == 0
                    kpathh=generateAlphabetCellArrayy(Knumber);
                end
                varargout=cell(1,6);
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
                varargout{6}=kcord;
            case 'automatic'
                kpathh=generateAlphabetCellArrayy(2);
                varargout=cell(1,6);
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
                varargout{6}=kcord;
            case 'crystal'
                if exist('kpathh','var') == 0
                    kpathh=generateAlphabetCellArrayy(2);
                end
                varargout=cell(1,6);
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
                varargout{6}=kcord;
            case 'tpiba'
                if exist('kpathh','var') == 0
                    kpathh=generateAlphabetCellArrayy(2);
                end
                varargout=cell(1,6);
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
                varargout{6}=kcord;
            case 'gamma'
                varargout=cell(1,6);
                kpathh{1}='Γ';
                varargout{1}=kpathh; 
                varargout{2}=kpathtickk_org;
                varargout{3}=kpathtickk_array;
                varargout{4}=kpathtickk;
                varargout{5}=kp_pointss;
                varargout{6}=kcord;
        end  
    end
    
    function cellArray = generateAlphabetCellArrayy(N)
        % Generate English letters
        letters = 'A':'Z';
        
        % Initialize cell array
        cellArray = cell(1,N);
        
        % Generate alphabet combinations
        index = 1;
        for i = 1:N
            % Calculate the index for letters array
            letterIndex = mod(index - 1, 26) + 1;
            
            % Calculate the number of repetitions for the prefix
            numRepetitions = ceil(index / 26);
            
            % Generate the prefix
            prefixx = repmat(letters(letterIndex), 1, numRepetitions);
            
            % Store the combination in the cell array
            cellArray{i} = prefixx;
            
            % Increment the index
            index = index + 1;
        end
    end

    function config = read_input_file(filename)
        % Default values
        default.WRITE_PDOS.prefix = [];
        default.WRITE_PDOS.inputdir = [];
        default.WRITE_PDOS.outdir = [];
        default.WRITE_PDOS.Nspin_Type = [];
        default.WRITE_PDOS.Calprefix = [];
        default.WRITE_PDOS.Efermi = [];
        default.WRITE_PDOS.Overwrite = [];
        default.WRITE_PDOS.Delimiter=[];
        
        default.PLOT_PDOS.Want_States = [];  % Cell array default
        default.PLOT_PDOS.Band_Limit = [];
        default.PLOT_PDOS.Color_Matrix = [];  % Default matrix
        default.PLOT_PDOS.Sorting_Type = [];  % Cell array default
        default.PLOT_PDOS.xy_limit = [];
        default.PLOT_PDOS.Interpoint = [];  % Default matrix
        default.PLOT_PDOS.Weight_Multiplier = [];  % Cell array default
        default.PLOT_PDOS.Justify_Zero = [];
        default.PLOT_PDOS.Figure_Spec = [];  % Default matrix
        default.PLOT_PDOS.Draw_Fermi = [];  % Cell array default
        default.PLOT_PDOS.Legend_Details = [];
        default.PLOT_PDOS.Aspect_Ratio = [];  % Default matrix
        default.PLOT_PDOS.Export_Figure = [];  % Default matrix
        default.PLOT_PDOS.Delimiter=[];
        default.PLOT_PDOS.Figure_Caption=[];
        % Initialize the config with default values
        config = default;
        
        % Open the file
        fid = fopen(filename, 'r');
        if fid == -1
            error('File could not be opened.');
        end
        
        current_section = '';
        
        % Read the file line by line
        while ~feof(fid)
            line = strtrim(fgetl(fid));
            % Remove comments: ignore anything after % or !
            comment_index = regexp(line, '[%!]', 'once');
            if ~isempty(comment_index)
                line = strtrim(line(1:comment_index-1));  % Trim line at comment
            end
            
            % Check for empty line after trimming
            if isempty(line)
                continue;
            end
            % Identify the section (&WRITE_PDOS, &PLOT_PDOS, etc.)
            if startsWith(line, '&')
                current_section = strrep(line, '&', '');  % Get section name
                continue;
            end
            
            % Identify the end of the section
            if strcmp(line, '/')
                current_section = '';
                continue;
            end
            
            % If we are in a section, process the key-value pairs
            if ~isempty(current_section)
                % Handle matrix variables like A=[1 2 3; 4 5 6],
                matrix_pattern = '(\w+)\s*=\s*\[(.+?)\],';
                match_matrix = regexp(line, matrix_pattern, 'tokens');
                if ~isempty(match_matrix)
                    key = match_matrix{1}{1};
                    value_str = match_matrix{1}{2};
                    
                    % Convert the matrix string into a MATLAB matrix
                    value = eval(['[' value_str ']']);

                    % Update config structure with the parsed value
                    config.(current_section).(key) = value;
                    continue;
                end
                
                % Handle cell array values like figure_spec={'Times','bold',12}
                cell_array_pattern = '(\w+)\s*=\s*\{(.+?)\},';
                match_cell = regexp(line, cell_array_pattern, 'tokens');
                if ~isempty(match_cell)
                    key = match_cell{1}{1};
                    value_str = match_cell{1}{2};
                    
                    % Convert the cell array string into a MATLAB cell array
                    value = eval(['{' value_str '}']);
                    
                    % Update config structure with the parsed value
                    config.(current_section).(key) = value;
                    continue;
                end
                
                % Handle normal key-value pairs like prefix="Niloy",
                tokens = regexp(line, '(\w+)\s*=\s*["'']?(.+?)["'']?,', 'tokens');
                if ~isempty(tokens)
                    key = tokens{1}{1};
                    value_str = tokens{1}{2};

                    % Attempt to convert the value to numeric if it's not a string
                    value_numeric = str2double(value_str);
                    if ~isnan(value_numeric)
                        value = value_numeric;  % Use numeric value if conversion is successful
                    else
                        value = value_str;  % Otherwise, keep it as a string
                    end
                    
                    % Update config structure with the parsed value
                    config.(current_section).(key) = value;
                end
            end
        end
        fclose(fid);
        
    end
end

