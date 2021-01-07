function sim_meta_fig_mofit(outputfolder,filename)
    
    csvfile=importdata(cat(2,outputfolder,filename));
    xmeans=csvfile(:,2);
    ymeans=csvfile(:,3);
    xvar=csvfile(:,4);
    yvar=csvfile(:,5);
    covar=csvfile(:,6);

    hamburger=figure('visible','off','Position', [10 10 1200 900]);
    
    subplot(1,3,1)
    title('Motility: Bivariate Normal Distribution Mean')
    boxplot([xmeans,ymeans],'Labels',{'X','Y'})  
    ylabel('Mean values of displacement per day.')
    
    subplot(1,3,2)
    title('Motility: Bivariate Normal Distribution Variance')
    boxplot([xvar,yvar],'Labels',{'X','Y'});
    ylabel('Variance values of displacement per day.')
    
    subplot(1,3,3)
    title('Motility: Bivariate Normal Distribution Covariance')
    boxplot(covar,'Labels',{'covariance'});
    ylabel('Covariance values of displacement per day.')
    
    saveas(hamburger,strcat(outputfolder,'analysis_meta_motility_fit_parameters.png'));
    close(hamburger)
    
end