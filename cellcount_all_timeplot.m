function maxtime=cellcount_all_timeplot(outputfolder,folder,timefile,identifier)
    % This function plots all cell populations.

    timemat=dlmread(strcat(folder,timefile), '\t', 1, 0)';
    fprintf(1,'Making full population plot.\n')
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    plot(timemat(1,:),timemat(7,:)+timemat(9,:)+timemat(11,:),...
         timemat(1,:),timemat(19,:)+timemat(21,:)+timemat(23,:),...
         timemat(1,:),timemat(13,:)+timemat(15,:)+timemat(17,:),...
         timemat(1,:),timemat(25,:),...
         timemat(1,:),timemat(7,:)+timemat(9,:)+timemat(11,:)+...
                      timemat(19,:)+timemat(21,:)+timemat(23,:)+...
                      timemat(13,:)+timemat(15,:)+timemat(17,:)+...
                      timemat(25,:)...
        );
    legend({'APC','CD8','CD4','Reticular','total'},'Location','north')
    title('Cell populations over time')
    xlabel("Time");
    ylabel("Cell number");
    saveas(humbolt,strcat(outputfolder,'cellcount_all_timeplot_',num2str(identifier),'.png'));
    close(humbolt)
    maxtime=timemat(1,end);
    
end