function cellcount_cd4_timeplot(outputfolder,folder,timefile,identifier)
     % This function plots all CD4 subpopulations.

     timemat=dlmread(strcat(folder,timefile), '\t', 1, 0)';
     fprintf(1,'Making CD4 population plot.\n')
     humbolt=figure('visible','off','Position', [10 10 1200 900]);
     plot(timemat(1,:),timemat(13,:)+timemat(15,:)+timemat(17,:),...
                timemat(1,:),timemat(13,:),...
                timemat(1,:),timemat(15,:),...
                timemat(1,:),timemat(17,:)...
         );
     legend({'total','Activated','Cognate','Non-cognate'},'Location','north')
     title('CD4 populations over time')
     xlabel("Time");
     ylabel("Cell number");
     saveas(humbolt,strcat(outputfolder,'cellcount_cd4_timeplot_',num2str(identifier),'.png'));
     close(humbolt)
    
%      dodo=figure('visible','off','Position', [10 10 1200 900]);
%      plot(timemat(1,:),timemat(13,:),...
%                 timemat(1,:),timemat(15,:),...
%                 timemat(1,:),timemat(13,:)+timemat(15,:)...
%        );
%      legend({'Activated','Cognate','Cognate total'},'Location','north')
%      title('Cognate CD4 population over time')
%      xlabel("Time");
%      ylabel("Cell number");
%      saveas(dodo,strcat(outputfolder,'cellcount_cd4ca_timeplot_',num2str(identifier),'.png'));
%      close(dodo)
    
end