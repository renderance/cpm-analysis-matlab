function motility_alt_measure_displacement(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)

    % Determines x and y displacement of each cell for each time-point.
    % Goes over all the ids, reads location data out of file and writes matrices for both coordinates.
    % It then calculates the displacement.
    
    
    fprintf(1,'Plotting alternative measure of displacement.\n');

    iterant=1:length(cellfiles)-1;
    time=(1:maxtime/10);
    xdelta=bsxfun(@minus,xcos(time(61:end),iterant),xcos(time(1:end-60),iterant));
    ydelta=bsxfun(@minus,ycos(time(61:end),iterant),ycos(time(1:end-60),iterant));
    zdelta=sqrt(bsxfun(@plus,ydelta.^2,xdelta.^2));
    vals=zdelta./(2*sqrt(750/pi));
    
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    boxplot(vals(1:(4*60):end,:)','BoxStyle','filled','MedianStyle','line','OutlierSize',2,'Symbol','.','labels',61:(4*60):maxtime/10);
    title('$\sqrt{(X_{n}-X_{n-1})^{2}+(Y_{n}-Y_{n-1})^{2}}$ displacement measure', 'Interpreter','latex');
    xlabel('Simulation Time (in minutes)')
    xtickangle(45)
    ylabel('Distance travelled per hour (in Cell diameters)') 
    saveas(humbolt,strcat(outputfolder,'motility_alt_measure_displacement_',num2str(identifier),'.png'));
    close(humbolt)
    
end