function barcodemat=barcode_time_extractor(barcodemat,folder,timefile)
    
    fprintf(1,'Working to extract interesting simulation properties.\n');
    form='%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n';
    factsarray=[];
    
    for ii=1:length(barcodemat(:,2))
        
        factarray=[];
        fileID=fopen(strcat(folder,timefile));
        
        if fileID == -1
        
            error('Called to non-existent file. Exiting.')
        
        else
            
            fgetl(fileID);
            line=fgetl(fileID);
            info=textscan(line,form);
            factarray=[info{4} info{5} info{20} info{24}]; 
            
        end
        
        factsarray=cat(1,factsarray,factarray);
        
    end
    
    barcodemat=cat(2,barcodemat,factsarray);

end