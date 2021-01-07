function motility_all_traceplot(outputfolder,xcos,ycos,famv,identifier)
    % This function creates a cell trace with unified start points of the last 2000 timesteps.

    % Determines x and y coordinates of each cell for each time-point within the last 2000.
    % Goes over all the ids, reads location data out of file and writes matrices for both coordinates.
    
    fprintf(1,'Making traceplot for all cells.\n');

    %pick=randi([1 length(cellfiles)-1],1,100);
    xtrace=xcos(:,:)-xcos(1,famv(:));
    ytrace=ycos(:,:)-ycos(1,famv(:));
    xtrace=xtrace./(2*sqrt(750/pi));
    ytrace=ytrace./(2*sqrt(750/pi));

    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    plot(xtrace(:,:),ytrace(:,:));
    title('Spatial trace for all cells');
    xlabel('Distance in X-direction (in cell diameters)');
    ylabel('Distance in Y-direction (in cell diameters)');
    saveas(humbolt,strcat(outputfolder,'motility_all_traceplot_',num2str(identifier),'.png'));
    close(humbolt)
    
end