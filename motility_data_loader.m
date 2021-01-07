function [xcos,ycos,famv,idsv,momv,barc,terv,actv,conv,motv]=motility_data_loader(cellfiles,maxtime)

    fprintf(1,'Loading motility data.\n')

    xcos=NaN(maxtime/10+1,length(cellfiles));
    ycos=NaN(maxtime/10+1,length(cellfiles));
    
    idsv=NaN(length(cellfiles),1);
    famv=NaN(length(cellfiles),1);
    momv=NaN(length(cellfiles),1);
    
    barc=[];
    
    terv=NaN(length(cellfiles),1);
    actv=NaN(length(cellfiles),1);
    conv=NaN(length(cellfiles),1);
    
    thing=-1;
    
    for iterant=1:length(cellfiles) 

        percentagefin=iterant/(length(cellfiles))*100;
        
        if floor(rem(percentagefin,10))==0 && thing<floor(percentagefin/10)
            
            fprintf(1,'\tFinished %d %% - proceeding with file %d of %d\n',int32(percentagefin),iterant,length(cellfiles));
            thing=thing+1;
            
        end
            
        filename=cellfiles(iterant);
        data=NaN(maxtime/10+1,4);
        data(:,1)=0:10:maxtime;
        temp=importdata(filename{1})';
        data(temp.data(1,1)/10+1:temp.data(end,1)/10+1,[2 3])=temp.data(:,[18 19]);
        xcos(:,temp.data(1,2))=data(:,2);
        ycos(:,temp.data(1,2))=data(:,3);
        
        idsv(temp.data(1,2))=temp.data(1,2);
        famv(temp.data(1,2))=temp.data(1,4);
        momv(temp.data(1,2))=temp.data(1,5);
        
        momv(idsv(and(famv~=idsv,and(famv~=0,momv==0))))=famv(idsv(and(famv~=idsv,and(famv~=0,momv==0))));
        %Ugly hack. Assumes that all cells from the family that do not have a known parent automatically hail from the naive cell.
        
        if ~ismember(temp.data(1,4),barc)
            
            barc=cat(1,barc,temp.data(1,4));
            
        end
        
        if temp.data(end,1)==maxtime
            
            terv(temp.data(1,2))=1;
            
        else
            
            terv(temp.data(1,2))=0;
            
        end
        
        if temp.data(end,3)==6
            
            actv(temp.data(1,2))=1;
            
        else
            
            actv(temp.data(1,2))=0;
            
        end
        
        prevstat=0;
        
        for zz=length(temp.data(:,13))
            
            newstat=temp.data(zz,13);
            
            if prevstat==0 && newstat~=0
                
                conv(temp.data(1,2))=conv(temp.data(1,2))+1;
                
            end
            
            prevstat=newstat;
            
        end
        
        clear data
        clear temp
        
    end
    
    motv=momv;
    motv(momv==0)=idsv(momv==0);
    
end