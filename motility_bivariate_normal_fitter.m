function [xmn,ymn,xvr,yvr,cvr]=motility_bivariate_normal_fitter(outputfolder,xcos,ycos,maxtime,identifier)

    fprintf(1,'Fitting travel distance for every day to bivariate normal distribution.\n');
    
    ydat=bsxfun(@minus,ycos(1441:1440:(maxtime/10),:),ycos(1:1440:(maxtime/10)-1440,:));
    xdat=bsxfun(@minus,xcos(1441:1440:(maxtime/10),:),xcos(1:1440:(maxtime/10)-1440,:));
    sise=size(ydat(:));
    ydis=reshape(ydat,sise);
    xdis=reshape(xdat,sise);
    ydis=ydis(~isnan(ydis));
    xdis=xdis(~isnan(xdis));
    
    humbolt=figure('visible','off','Position', [10 10 1200 900]);
    scatter(ydis./(2*sqrt(750/pi)),xdis./(2*sqrt(750/pi)),2,'.','MarkerFaceColor','k');
    
    xydis=cat(2,xdis,ydis);
    xydis=xydis./(2*sqrt(750/pi));
    
    xymn=nanmean(xydis);
    xycv=cov(xydis);
    
    xsteps = linspace(min(xydis(:,1)),max(xydis(:,1)));
    ysteps = linspace(min(xydis(:,2)),max(xydis(:,2)));
    [X,Y] = meshgrid(xsteps,ysteps);
    F = mvnpdf([X(:) Y(:)],xymn,xycv);
    
    F = reshape(F,length(ysteps),length(xsteps));
    hold on
    contour(xsteps,ysteps,F);
    
    title('Cell displacement in a day.');
    xlabel('Displacement in X-direction (in cell diameter)')
    ylabel('Displacement in Y-direction (in cell diameter)')
    saveas(humbolt,strcat(outputfolder,'motility_bivariate_normal_fit_',num2str(identifier),'.png'));
    close(humbolt)
    
    xmn=xymn(1);
    ymn=xymn(2);
    xvr=xycv(1,1);
    yvr=xycv(2,2);
    cvr=xycv(1,2);
    
end