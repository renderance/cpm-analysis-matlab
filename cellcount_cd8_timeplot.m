function cellcount_cd8_timeplot(outputfolder,folder,timefile,identifier)
    % This function plots CD8 subpopulations.

    timemat=dlmread(strcat(folder,timefile), '\t', 1, 0)';
    fprintf(1,'Making CD8 population plot.\n')
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    plot(timemat(1,:),timemat(19,:)+timemat(21,:)+timemat(23,:),...
               timemat(1,:),timemat(19,:),...
               timemat(1,:),timemat(21,:),...
               timemat(1,:),timemat(23,:)...
        );
    legend({'total','Activated','Cognate','Non-cognate'},'Location','north')
    title('CD8 populations over time')
    xlabel("Time");
    ylabel("Cell number");
    saveas(humbolt,strcat(outputfolder,'cellcount_cd8_timeplot_',num2str(identifier),'.png'));
    close(humbolt)

    dodo=figure('visible','off','Position', [10 10 1200 900]);
    plot(timemat(1,:),timemat(19,:),...
               timemat(1,:),timemat(21,:),...
               timemat(1,:),timemat(19,:)+timemat(21,:)...
       );
    legend({'Activated','Cognate','Cognate total'},'Location','north')
    title('Cognate CD8 population over time')
    xlabel("Time");
    ylabel("Cell number");
    saveas(dodo,strcat(outputfolder,'cellcount_cd8ca_timeplot_',num2str(identifier),'.png'));
    close(dodo)

end