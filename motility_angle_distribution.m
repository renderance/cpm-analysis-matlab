function motility_angle_distribution(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)

    fprintf(1,'Making polar histogram for all cell angles.\n');
    
    iterant=1:length(cellfiles)-1;
    time=(1:maxtime/10);
    xdelta=bsxfun(@minus,xcos(time+1,iterant),xcos(time,iterant));
    ydelta=bsxfun(@minus,ycos(time+1,iterant),ycos(time,iterant));
    
    anglesr=atan2(xdelta,ydelta)+pi;
    anglesd=radtodeg(anglesr-pi)+180;
    
    counts=NaN(1,360);
    
    counts(1)=sum(anglesd(:)<=0.5);
    for ii=2:360
        counts(ii)=sum(anglesd(:)<=ii-0.5)-sum(counts(1:ii-1));
    end
    counts(1)=counts(1)+sum(anglesd(:)>359.5);
    
    binedges=-0.5:1:360-0.5;
    polarbinedges=binedges.*(2*pi/360);
    
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    polarhistogram('BinEdges',polarbinedges,'BinCounts',counts);
    thetaticks(0:45:315);
    title('Cell directions as angles relative to field x-axis');
    saveas(humbolt,strcat(outputfolder,'motility_angle_distribution_polar_',num2str(identifier),'.png'));
    close(humbolt)
    
    dodo=figure('visible','off','Position', [10 10 1200 900]);
    turf=histogram('BinEdges',binedges,'BinCounts',counts);
    set(turf.Parent,'YScale','log');
    xticks(0:45:360);
    title('Cell directions as angles relative to field x-axis');
    saveas(dodo,strcat(outputfolder,'motility_angle_distribution_',num2str(identifier),'.png'));
    close(dodo)

end
