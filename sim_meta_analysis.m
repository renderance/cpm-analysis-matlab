function sim_meta_analysis(parentfolder,outputfolder,targettime,last)

    [status,string]=system(cat(2,'ls -1 ',char(parentfolder),'sim_*.tar.gz'));
    
    if status~=0
        
        fprintf('Error: Could not list simulations in this parent folder.\n')
        
    end
    
    newstring=splitlines(string);
    newstring=newstring(1:end-1);
    
    tafolder='temp_analysis_folder';
    
    if last == 1
        
        system(cat(2,'rm -r ',char(outputfolder),tafolder));
        system(cat(2,'rm -r ',char(outputfolder),'*.csv'));
        system(cat(2,'rm -r ',char(outputfolder),'*.png'));
        
    end
    
    for ii=last:length(newstring)
        
        folder=newstring{ii};
        identifier=str2double(folder(length(parentfolder)+1+length('sim_'):end-length('.tar.gz')));
        
        if isnan(identifier)
            
            error('Error: Simulation with NaN identifier detected.\n');
            
        end
        
        fprintf('Working on analyis and data-extraction from simulation %d.\n',identifier);
        
        stat1=system(cat(2,'mkdir ',char(outputfolder),tafolder));
        
        if stat1~=0
           
            fprintf('Error: Could not make analysis folder.\n');
            
        end
        
        stat2=system(cat(2,'tar -xf ',folder,' -C ',char(outputfolder),tafolder));
        
        if stat2~=0
            
            fprintf('Error: Could not untar.\n');
            
        end
        
        fulltafolder=cat(2,char(outputfolder),tafolder);
        
        [status2,string2]=system(cat(2,'ls -1 ',char(fulltafolder),'/cell_*.csv'));
        
        if status2~=0
            
            fprintf('Error: Could not list cell files.\n');
            
        end
        
        newstring2=splitlines(string2);
        cellfiles=newstring2(1:end-1);
        timefile='/time.csv';
               
        [xmean,ymean,xvar,yvar,covar,fstat,undiv,div,all,bcmat]=...
            sim_analysis_combiner(outputfolder,fulltafolder,cellfiles,timefile,identifier,targettime);
        
        bcmat=cat(2,bcmat(:,[1 2 3]),size(bcmat,1).*ones(size(bcmat,1),1),bcmat(:,4:end));
        
        line1=[identifier,xmean,ymean,xvar,yvar,covar];
        dlmwrite(cat(2,outputfolder,'analysis_meta_sim_motility_data.csv'),line1,'delimiter','\t','-append');
        
        line2=[identifier,fstat,undiv,div,all];
        dlmwrite(cat(2,outputfolder,'analysis_meta_sim_division_status.csv'),line2,'delimiter','\t','-append');
        
        dlmwrite(cat(2,outputfolder,'analysis_meta_barcode_data.csv'),bcmat,'delimiter','\t','-append');
        
        stat3=system(cat(2,'rm -r ',char(outputfolder),tafolder));
        
        if stat3~=0
            
            fprintf('Error: Could not remove analysis folder.\n');
            
        end
        
        if ii > 3
            
            sim_meta_fig_mofit(outputfolder,'analysis_meta_sim_motility_data.csv');
            sim_meta_fig_divfrac(outputfolder,'analysis_meta_sim_division_status.csv');
            % meta_analysis_figmaker;
            
        end
        
        fclose('all');
                
    end
    
end