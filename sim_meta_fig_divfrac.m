function sim_meta_fig_divfrac(outputfolder,filename)
    
    csvfile=importdata(cat(2,outputfolder,filename));
    undvec=csvfile(:,3);
    divvec=csvfile(:,4);
    allvec=csvfile(:,5);

    hamburger=figure('visible','off','Position', [10 10 1200 900]);
    
    title('Division behaviour.')
    %boxplot([undvec./allvec,divvec./allvec],'labels',{'Undivided','Divided'})  
    histogram(undvec./allvec,-0.05:0.1:1.05)
    xlabel('Fraction of cells that divided.')
    ylim([0 size(undvec,1)]);
    ylabel('Number of simulations')
    
    saveas(hamburger,strcat(outputfolder,'analysis_meta_division_fates.png'));
    close(hamburger)
    
end