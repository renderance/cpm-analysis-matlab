function motility_turning_angle_distribution(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)

    fprintf(1,'Making polar histogram for all turning angles.\n');
    
    iterant=1:length(cellfiles)-1;
    time=(1:maxtime/10);
    xdelta=bsxfun(@minus,xcos(time(16:end),iterant),xcos(time(1:end-15),iterant));
    ydelta=bsxfun(@minus,ycos(time(16:end),iterant),ycos(time(1:end-15),iterant));
    angles=atan2(xdelta,ydelta);
    
    radangx=angles(1:end-1,:);
    radangy=angles(2:end,:);
    radangxv=radangx(:);
    radangyv=radangy(:);
    
    cosradangxv=cos(radangxv);
    cosradangyv=cos(radangyv);
    sinradangxv=sin(radangxv);
    sinradangyv=sin(radangyv);
    
    partcostur=cosradangxv.*cosradangyv;
    partsintur=sinradangxv.*sinradangyv;
    dubtur=bsxfun(@plus,partcostur,partsintur);
    turnur=real(acos(dubtur));
    
    turnan=radtodeg(turnur);
    %turnan(turnan<0)=360+turnan(turnan(turnan<0));
    
    turning=reshape(turnan,length(time(1:end-16)),length(iterant));
    
    counts=NaN(1,360);
    
    counts(1)=sum(turning(:)<=1);
    for ii=2:360
        counts(ii)=sum(turning(:)<=ii)-sum(counts(1:ii-1));
    end
    counts(1)=counts(1)+sum(turning(:)>360);
    
    binedges=0:1:360;
    polarbinedges=binedges.*(2*pi/360);
    
    
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    polarhistogram('BinEdges',polarbinedges,'BinCounts',counts);
    title('15-minute turning angle distribution');
    saveas(humbolt,strcat(outputfolder,'motility_turning_angle_distribution_polar_',num2str(identifier),'.png'));
    close(humbolt)
    
    dodo=figure('visible','off','Position', [10 10 1200 900]);
    turf=histogram('BinEdges',0:1:180,'BinCounts',counts(1:180));
    set(turf.Parent,'YScale','log');
    title('15-minute turning angle distribution');
    saveas(dodo,strcat(outputfolder,'motility_turning_angle_distribution_',num2str(identifier),'.png'));
    close(dodo)

end