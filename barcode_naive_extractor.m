function barcodes=barcode_naive_extractor(barcodemat,folder,maxtime)

    fprintf(1,'Working to extract interesting naive properties.\n');
    form='%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n';
    factsarray=[];
    
    for ii=1:length(barcodemat(:,2))
        
        factarray=[];
        fileID=fopen(strcat(folder,'/cell_',num2str(barcodemat(ii,2)),'.csv'));
        
        if fileID == -1
        
            error('Called to non-existent file. Exiting.')
        
        else
            
            fgetl(fileID);
            line=fgetl(fileID);
            info=textscan(line,form);
            factarray=cat(2,factarray,[info{16} info{17} info{8} info{14} info{9}]);
            
            activated=0;
            broke=0;
            left=0;
            prevstat=0;
            actevent=0;
            
            actinfo=    [NaN NaN NaN NaN];
            brinfo=     [NaN NaN NaN NaN];
            leavinfo=   [NaN NaN NaN NaN];
            divinfo=    [NaN NaN NaN NaN];
            
            
            while feof(fileID)~=1
                    
                line=fgetl(fileID);
                info=textscan(line,form);
            
                if info{3}~=7 && activated==0
                
                    actinfo=[info{1} info{10} info{12} info{15}];
                    activated=1;
                    
                end
                
                if info{11}==0 && broke==0
                    
                    
                    brinfo=[info{1} info{10} info{12} info{15}];
                    broke=1;
                    
                end
                
                if info{3}~=7 && info{11}==0 && left==0
                    
                    leavinfo=[info{1}-actinfo(1) info{10} info{12} info{15}];
                    left=1;
                    
                end
                
                newstat=info{11};
                
                if prevstat==0 && newstat==1
                    
                    actevent=actevent+1;
                    
                end
                
                prevstat=newstat;
                
            end
            
            info=textscan(line,form);
            
            if info{1}~=maxtime
                
                divinfo=[info{1} info{10} info{12} info{15}];
                
            end
            
            factarray=cat(2,factarray,actinfo,brinfo,divinfo,leavinfo,actevent);
              
        end
            
        factsarray=cat(1,factsarray,factarray);
        fclose(fileID);
        
    end
    
    barcodes=cat(2,barcodemat,factsarray);

end