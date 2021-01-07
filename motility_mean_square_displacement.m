function motility_mean_square_displacement(outputfolder,xcos,ycos,famv,cellfiles,maxtime,identifier)

    % Determines x and y displacement of each cell for each time-point.
    % Goes over all the ids, reads location data out of file and writes matrices for both coordinates.
    % It then calculates the displacement.
    
    fprintf(1,'Plotting mean square displacement.\n');

    xtrace=bsxfun(@minus,xcos(:,:),xcos(1,famv(:)));
    ytrace=bsxfun(@minus,ycos(:,:),ycos(1,famv(:)));
    ztrace=sqrt(bsxfun(@plus,ytrace.^2,xtrace.^2));
    
    what=31:30:(maxtime/10)+1;
    msds=NaN(30,length(cellfiles)-1);
    
    for jterant=1:30
        
        for kterant=1:length(cellfiles)-1
        
            msds(jterant,kterant)=...
                1/length(what)*...
                sum((ztrace(what,kterant)-ztrace(what-jterant,kterant)).^2);
            
        end
        
    end
    
    vals=msds/750;
    
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    boxplot(vals(:,:)','BoxStyle','filled','MedianStyle','line','OutlierSize',2,'Symbol','.','labels',1:30);
    title('Mean square displacement');
    xlabel('Time (Minutes)');
    ylabel('Mean Square Displacement (in T-cell volumes)');
    
    % Mean square displacement of each cell, calculated from length(what) 
    % number of samples. Each sample was taken from each cell's individual
    % trajectory, and from each sample, sub-trajectories of the used 
    % lenghts were taken to calculate MSD's. The mean was thus taken over 
    % multiple instances of the same cell's total simulated migration, and
    % calculated for increasingly long displacement periods.
    
    saveas(humbolt,strcat(outputfolder,'motility_mean_square_displacement_',num2str(identifier),'.png'));
    close(humbolt)
    
end