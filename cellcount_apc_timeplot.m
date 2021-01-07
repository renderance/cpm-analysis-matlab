function cellcount_apc_timeplot(outputfolder,folder,timefile,identifier)
    % This function plots all APC subpopulations.

    timemat=dlmread(strcat(folder,timefile), '\t', 1, 0)';
    fprintf(1,'Making APC population plot.\n')
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    plot(timemat(1,:),timemat(7,:)+timemat(9,:)+timemat(11,:),...
               timemat(1,:),timemat(7,:),...
               timemat(1,:),timemat(9,:),...
               timemat(1,:),timemat(11,:)...
        );
    legend({'total','Activated','Cognate','Non-cognate'},'Location','north')
    title('APC populations over time')
    xlabel("Time");
    ylabel("Cell number");
    saveas(humbolt,strcat(outputfolder,'cellcount_apc_timeplot_',num2str(identifier),'.png'));
    close(humbolt)
    
end