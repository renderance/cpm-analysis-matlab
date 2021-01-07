
stat0=system(cat(2,'mkdir ','/home/glenn/Documents/morpheus/sim_outputs3'));
stat1=system(cat(2,'mkdir ','/home/glenn/Documents/morpheus/sim_outputs3/temp_analysis_folder'));
stat2=system(cat(2,'tar -xf ','/hosts/linuxhome/vacuole2/tmp/glenn2/sim_51030011.tar.gz',' -C ','/home/glenn/Documents/morpheus/sim_outputs3/','temp_analysis_folder'));
set(0,'DefaultFigureWindowStyle','docked')
if stat2~=0
    fprintf('Error: Could not untar.\n');
end
fulltafolder=cat(2,'/home/glenn/Documents/morpheus/sim_outputs3/','temp_analysis_folder');
[status2,string2]=system(cat(2,'ls -1 ',char(fulltafolder),'/cell_*.csv'));
if status2~=0
    fprintf('Error: Could not list cell files.\n');
end
newstring2=splitlines(string2);
cellfiles=newstring2(1:end-1);
timefile='/time.csv';
maxtime=72000;
[xcos,ycos,famv,idsv,momv,~,~,~,~,motv]=motility_data_loader(cellfiles,maxtime);

time=1:7200;
iterant=1:length(cellfiles)-1;
xdelta=bsxfun(@minus,xcos(time+1,iterant),xcos(time,iterant));
ydelta=bsxfun(@minus,ycos(time+1,iterant),ycos(time,iterant));
anglesr=atan2(xdelta,ydelta)+pi;
anglesd=radtodeg(anglesr-pi)+180;

mat=zeros(361,720);
for jterant=1:7199
    for zterant=1:length(cellfiles)-1
        anglet=anglesd(jterant,zterant);
        mat(round(anglet,0)+1,floor(jterant/10)+1)=mat(round(anglet,0)+1,floor(jterant/10)+1)+1;
    end
end
mat(1,:)=mat(1,:)+mat(361,:);
mat=mat(1:360,:);
for vterant=1:720
    mat(:,vterant)=mat(:,vterant)./sum(mat(:,vterant));
end













figure; % ORIENTATION PLOT
imagesc(mat);
xticks(0:144:720);
%xtickangle(270);
xticklabels(0:1:5)
xlabel('Time (in days)')
yticks([0.5:45:315.5 359.5])
yticklabels([0:45:315 359])
ylabel('Orientation (in degrees)')
colorbar;
title("Cell orientation");
caxis([0 0.01]);

tanglesd2=bsxfun(@minus,anglesd(time(2:end),iterant),anglesd(time(1:end-1),iterant));
tanglesd4=bsxfun(@minus,anglesd(time(4:end),iterant),anglesd(time(1:end-3),iterant));
tanglesd8=bsxfun(@minus,anglesd(time(8:end),iterant),anglesd(time(1:end-7),iterant));
tanglesd16=bsxfun(@minus,anglesd(time(16:end),iterant),anglesd(time(1:end-15),iterant));
tanglesd32=bsxfun(@minus,anglesd(time(32:end),iterant),anglesd(time(1:end-31),iterant));
tanglesd64=bsxfun(@minus,anglesd(time(64:end),iterant),anglesd(time(1:end-63),iterant));
tanglesd128=bsxfun(@minus,anglesd(time(128:end),iterant),anglesd(time(1:end-127),iterant));
tanglesd256=bsxfun(@minus,anglesd(time(256:end),iterant),anglesd(time(1:end-255),iterant));

countsd2=zeros(360,1);
countsd4=zeros(360,1);
countsd8=zeros(360,1);
countsd16=zeros(360,1);
countsd32=zeros(360,1);
countsd64=zeros(360,1);
countsd128=zeros(360,1);
countsd256=zeros(360,1);

for zterant=1:size(tanglesd2,2)
    for qterant=1:360
        countsd2(qterant)=countsd2(qterant)+sum(round(tanglesd2(:,zterant))+1==qterant);
        countsd4(qterant)=countsd4(qterant)+sum(round(tanglesd4(:,zterant))+1==qterant);
        countsd8(qterant)=countsd8(qterant)+sum(round(tanglesd8(:,zterant))+1==qterant);
        countsd16(qterant)=countsd16(qterant)+sum(round(tanglesd16(:,zterant))+1==qterant);
        countsd32(qterant)=countsd32(qterant)+sum(round(tanglesd32(:,zterant))+1==qterant);
        countsd64(qterant)=countsd64(qterant)+sum(round(tanglesd64(:,zterant))+1==qterant);
        countsd128(qterant)=countsd128(qterant)+sum(round(tanglesd128(:,zterant))+1==qterant);
        countsd256(qterant)=countsd256(qterant)+sum(round(tanglesd256(:,zterant))+1==qterant);
        fprintf('Num:\t%d\t%d\n',zterant,qterant);
    end
    countsd2(1)=countsd2(1)+sum(round(tanglesd2(:,zterant))+1==361);
    countsd4(1)=countsd4(1)+sum(round(tanglesd4(:,zterant))+1==361);
    countsd8(1)=countsd8(1)+sum(round(tanglesd8(:,zterant))+1==361);
    countsd16(1)=countsd16(1)+sum(round(tanglesd16(:,zterant))+1==361);
    countsd32(1)=countsd32(1)+sum(round(tanglesd32(:,zterant))+1==361);
    countsd64(1)=countsd64(1)+sum(round(tanglesd64(:,zterant))+1==361);
    countsd128(1)=countsd128(1)+sum(round(tanglesd128(:,zterant))+1==361);
    countsd256(1)=countsd256(1)+sum(round(tanglesd256(:,zterant))+1==361);
    fprintf('Num:\t%d\t%d\n',zterant,qterant+1);
end
counts=cat(2,countsd2,countsd4,countsd8,countsd16,countsd32,countsd64,countsd128,countsd256);
for cterant=1:8
    counts(:,cterant)=counts(:,cterant)./sum(counts(:,cterant));
end

countsd2m=zeros(size(countsd2));
countsd4m=zeros(size(countsd4));
countsd8m=zeros(size(countsd8));
countsd16m=zeros(size(countsd16));
countsd32m=zeros(size(countsd32));
countsd64m=zeros(size(countsd64));
countsd128m=zeros(size(countsd128));
countsd256m=zeros(size(countsd256));

for wterant=1:180
    countsd2m(wterant)=countsd2(wterant)+countsd2(361-wterant);
    countsd4m(wterant)=countsd4(wterant)+countsd4(361-wterant);
    countsd8m(wterant)=countsd8(wterant)+countsd8(361-wterant);
    countsd16m(wterant)=countsd16(wterant)+countsd16(361-wterant);
    countsd32m(wterant)=countsd32(wterant)+countsd32(361-wterant);
    countsd64m(wterant)=countsd64(wterant)+countsd64(361-wterant);
    countsd128m(wterant)=countsd128(wterant)+countsd128(361-wterant);
    countsd256m(wterant)=countsd256(wterant)+countsd256(361-wterant);
end
counts=cat(2,countsd2m,countsd4m,countsd8m,countsd16m,countsd32m,countsd64m,countsd128m,countsd256m);
for cterant=1:8
    counts(:,cterant)=counts(:,cterant)./sum(counts(:,cterant));
end



















figure; % TURNING ANGLE
stairs(axes,counts)
xticks([0.5:22.5:180.5])
xticklabels([0:22.5:180])
set(gca,'YScale','log');
ylabel('Observations (fraction of total)');
xlabel('Angle (in degrees)');
title("Turning angle distribution");
legend('2 min','4 min','8 min','16 min','32 min','1 hour 4 min')
xlim([0.5 180.5])

xtrace=bsxfun(@minus,xcos(:,:),xcos(1,motv(:)));
ytrace=bsxfun(@minus,ycos(:,:),ycos(1,motv(:)));
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














figure; % MEAN SQUARE DISPLACEMENT
boxplot(vals(:,:)','BoxStyle','outline','MedianStyle','line','OutlierSize',2,'Symbol','.','labels',1:30);
title('Mean square displacement');
xlabel('Time (minutes)');
ylabel('Mean Square Displacement (in T-cell volumes)');
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')

ydat=bsxfun(@minus,ycos(1441:1440:(maxtime/10),:),ycos(1:1440:(maxtime/10)-1440,:));
xdat=bsxfun(@minus,xcos(1441:1440:(maxtime/10),:),xcos(1:1440:(maxtime/10)-1440,:));
sise=size(ydat(:));
ydis=reshape(ydat,sise);
xdis=reshape(xdat,sise);
ydis=ydis(~isnan(ydis));
xdis=xdis(~isnan(xdis));














figure; % BIVARIATE NORMAL FIT
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
colormap('hot')
contour(xsteps,ysteps,F);
title('Cell displacement in a day.');
xlabel('Displacement in X-direction (in cell diameter)')
ylabel('Displacement in Y-direction (in cell diameter)')

csvfile=importdata(cat(2,'/home/glenn/Documents/morpheus/sim_outputs2/','analysis_meta_sim_motility_data.csv'));
xmeans=csvfile(:,2);
ymeans=csvfile(:,3);
xvar=csvfile(:,4);
yvar=csvfile(:,5);
covar=csvfile(:,6);












figure; % VARIANCE COVARIANCE PLOT
subplot(1,3,1)
title('Motility: Bivariate Normal Distribution Mean')
boxplot([xmeans,ymeans],'Labels',{'X','Y'})  
ylabel('Mean values of displacement per day (in gridspaces).')
subplot(1,3,2)
title('Motility: Bivariate Normal Distribution Variance')
boxplot([xvar,yvar],'Labels',{'X','Y'});
ylabel('Variance values of displacement per day (in gridspaces^2).')
subplot(1,3,3)
title('Motility: Bivariate Normal Distribution Covariance')
boxplot(covar,'Labels',{'covariance'});
ylabel('Covariance values of displacement per day (in gridspaces^2).')
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')



















% IMPORT DATA FOR OTHER FIGURES
bcd1=importdata('/home/glenn/Documents/morpheus/sim_outputs/analysis_meta_barcode_data.csv');
bcs1=sortrows(bcd1,[3 33],'descend');
bcd2=importdata('/home/glenn/Documents/morpheus/sim_outputs2/analysis_meta_barcode_data.csv');
bcs2=sortrows(bcd2,[3 33],'descend');

cpv1=bcd1(:,28)./0.250;
cpv2=bcd2(:,28)./0.250;
cpv3=((1./bcd1(:,27)).*4000000)./0.250;
cpv4=((1./bcd2(:,27)).*4000000)./0.250;













f1=figure; % CELL DENSITY DISTRIBUTION
a1=axes;
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
h1=histogram(a1,cpv1,'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
hold on
h2=histogram(a1,cpv2,'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
xlim([10000 30000]);
xlabel('Cell density ($\frac{\mathrm{cells}}{\mathrm{mm}^2}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O');
title('Distribution of simulated starting cell densities','FontSize',34);













f1=figure; % SEEDED COGNATE CELLS
a1=axes;
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
h1=histogram(a1,bcd1(:,4),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
h1.Normalization='probability';
h1.FaceAlpha=0.4;
h1.FaceColor='b';
hold on
h2=histogram(a1,bcd2(:,4),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
h2.Normalization='probability';
xlim([0 15]);
h2.FaceAlpha=0.4;
h2.FaceColor='r';
xlabel('Number of cognate CD8$^{+}$ cells','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Starting cognate cell number','FontSize',34);














f1=figure; % CONTACT REQUIREMENT DISTRIBUTION
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,bcd1(:,7),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
hold on
h2=histogram(a1,bcd2(:,7),'DisplayName','Inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
xlabel('Contact requirement (Req$_{i}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Naive cell contact requirement','FontSize',34);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);













f1=figure; % STIMULATION STRENGTH
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,bcd1(:,9),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinEdges',0.6:0.025:1.4);
hold on
h2=histogram(a1,bcd2(:,9),'DisplayName','Inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinEdges',0.6:0.025:1.4);
xlabel('Dendritic cell effectiveness ($c_{\mathrm{D}}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Stimulation effectiveness','FontSize',34);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);













f1=figure; % SENSITIVITY PARAMETER DISTRIBUTION
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,bcd1(:,8),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
hold on
h2=histogram(a1,bcd2(:,8),'DisplayName','Inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
xlabel('Starting Myc (Myc$_{i,0}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Naive cell stimulation sensitivity','FontSize',34);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);















f1=figure; % FIRST DIVISION TIME
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,bcd1(:,5),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
hold on
h2=histogram(a1,bcd2(:,5),'DisplayName','Inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
xlabel('Time of first division ($t_{i,\mathrm{div},1}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex')
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Naive cell first division clock','FontSize',34);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);













f1=figure; % SUBSEQUENT DIVISION TIME
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,bcd1(:,6),'DisplayName','No inheritance','FaceAlpha',0.4,'FaceColor',[0 0 1],'Normalization','probability','BinMethod','auto');
hold on
h2=histogram(a1,bcd2(:,6),'DisplayName','Inheritance','FaceAlpha',0.4,'FaceColor',[1 0 0],'Normalization','probability','BinMethod','auto');
xlabel('Time of division ($t_{i,\mathrm{div}}$)','FontSize',26,'FontName','Linux Libertine Display O','Interpreter','latex');
ylabel('Probability','FontSize',32,'FontName','Linux Libertine Display O')
title('Naive cell division clock','FontSize',34);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);















% IGNORE UNDIVIDED CELLS FOR FUTURE PLOTS
dbcd1=bcd1(bcd1(:,3)>1,:);
dbcd2=bcd2(bcd2(:,3)>1,:);














f1=figure; % TIME TO FIRST DIVISION
edges=0:600:72000;
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,dbcd1(:,18),edges,'DisplayName','No inheritance','FaceAlpha',0.25,'FaceColor',[0 0 1],'Normalization','probability','EdgeAlpha',0.25);
xticks(0:14400:72000);
xticklabels(0:1:5);
xlabel(a1,'Time to first division (days)','FontSize',26,'FontName','Linux Libertine Display O')
xlim([0 72000])
ylabel(a1,'Frequency','FontSize',32,'FontName','Linux Libertine Display O')
hold on
h2=histogram(a1,dbcd2(:,18),edges,'DisplayName','Inheritance','FaceAlpha',0.25,'FaceColor',[1 0 0],'Normalization','probability','EdgeAlpha',0.25);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
title('Time to first division','FontSize',34,'FontName','Linux Libertine Display O');
%pd1=fitdist(dbcd1(:,18),'normal');
%pd2=fitdist(dbcd2(:,18),'normal');
%pdf1=pdf(pd1,dbcd1(:,18));
%pdf2=pdf(pd2,dbcd2(:,18));
%pdf1=pdf1*sum(h1.Values * h1.BinWidth);
%pdf2=pdf2*sum(h2.Values * h2.BinWidth);
%[x1, idx1]=sort(dbcd1(:,18));
%[x2, idx2]=sort(dbcd2(:,18));
%y1=pdf1(idx1);
%y2=pdf2(idx2);
%plot(x1,y1,'b-','LineWidth',1.5,'DisplayName','Fit: No Inheritance');
%plot(x2,y2,'r-','LineWidth',1.5,'DisplayName','Fit: Inheritance');













f1=figure; % ADDITIONAL ACTIVATION
edges=0:60:3600;
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,dbcd1(:,22),edges,'DisplayName','No inheritance','FaceAlpha',0.25,'FaceColor',[0 0 1],'Normalization','probability','EdgeAlpha',0.25);
xticks(0:600:3600);
xticklabels(0:1:6);
ylabel(a1,'Frequency','FontSize',32,'FontName','Linux Libertine Display O')
xlim([0 3600])
hold on
h2=histogram(a1,dbcd2(:,22),edges,'DisplayName','Inheritance','FaceAlpha',0.25,'FaceColor',[1 0 0],'Normalization','probability','EdgeAlpha',0.25);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
title('Additional activation','FontSize',34,'FontName','Linux Libertine Display O');
ylabel(a1,'Frequency','FontSize',32,'FontName','Linux Libertine Display O')
xlabel(a1,'Time between activation and break away from APC (hours)','FontSize',28,'FontName','Linux Libertine Display O');
frac1=sum(dbcd1(:,22)>1800)/length(dbcd1(:,22));
frac2=sum(dbcd2(:,22)>1800)/length(dbcd2(:,22));
%pd1=fitdist(dbcd1(:,22),'Weibull');
%pd2=fitdist(dbcd2(:,22),'Weibull');
%pdf1=pdf(pd1,dbcd1(:,22));
%pdf2=pdf(pd2,dbcd2(:,22));
%pdf1=pdf1*sum(h1.Values * h1.BinWidth);
%pdf2=pdf2*sum(h2.Values * h2.BinWidth);
%[x1, idx1]=sort(dbcd1(:,22));
%[x2, idx2]=sort(dbcd2(:,22));
%y1=pdf1(idx1);
%y2=pdf2(idx2);
%plot(x1,y1,'b-','LineWidth',1.5,'DisplayName','Fit: No Inheritance');
%plot(x2,y2,'r-','LineWidth',1.5,'DisplayName','Fit: Inheritance');












f1=figure; % TIME TO ACTIVATION
edges=0:60:6000;
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
h1=histogram(a1,dbcd1(:,10),edges,'DisplayName','No inheritance','FaceAlpha',0.25,'FaceColor',[0 0 1],'Normalization','probability','EdgeAlpha',0.25);
xticks(0:1200:6000);
xticklabels(0:2:10);
ylabel(a1,'Frequency','FontSize',32,'FontName','Linux Libertine Display O')
xlabel(a1,'Time to activation (hours)','FontSize',26,'FontName','Linux Libertine Display O')
xlim([0 6000])
hold on
h2=histogram(a1,dbcd2(:,10),edges,'DisplayName','Inheritance','FaceAlpha',0.25,'FaceColor',[1 0 0],'Normalization','probability','EdgeAlpha',0.25);
legend1 = legend(a1,'show');
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
set(a1,'FontName','Linux Libertine Display O','FontSize',22);
title('Time to activation','FontSize',34,'FontName','Linux Libertine Display O');
frac1=sum(dbcd1(:,10)>6000)/length(dbcd1(:,10));
frac2=sum(dbcd2(:,10)>6000)/length(dbcd2(:,10));
%pd1=fitdist(dbcd1(:,10),'loglogistic');
%pd2=fitdist(dbcd2(:,10),'loglogistic');
%pdf1=pdf(pd1,dbcd1(:,10));
%pdf2=pdf(pd2,dbcd2(:,10));
%pdf1=pdf1*sum(h1.Values * h1.BinWidth);
%pdf2=pdf2*sum(h2.Values * h2.BinWidth);
%[x1, idx1]=sort(dbcd1(:,10));
%[x2, idx2]=sort(dbcd2(:,10));
%y1=pdf1(idx1);
%y2=pdf2(idx2);
%plot(x1,y1,'b-','LineWidth',1.5,'DisplayName','Fit: No Inheritance');
%plot(x2,y2,'r-','LineWidth',1.5,'DisplayName','Fit: Inheritance');

























humphrey=figure; % BARCODE RANKING
timothy=subplot(1,2,1);
bar(timothy,bcs1(:,3),'b');
hold on
bar(timothy,bcs1(:,33),'r');
ylabel(timothy,'Family size / Number of active daughters','FontSize',22,'FontName','Linux Libertine Display O')
xlabel(timothy,'Barcodes','FontSize',22,'FontName','Linux Libertine Display O')
xtickangle(270)
xticks([]);
timothy.TickDir='out';
title('No Inheritance','FontSize',30,'FontName','Linux Libertine Display O');
set(timothy,'FontName','Linux Libertine Display O','FontSize',22,'yscale','log');
henry=subplot(1,2,2);
bar(henry,bcs2(:,3),'b');
hold on
bar(henry,bcs2(:,33),'r');
ylabel(henry,'Family size / Number of active daughters','FontSize',22,'FontName','Linux Libertine Display O')
xlabel(henry,'Barcodes','FontSize',22,'FontName','Linux Libertine Display O')
xtickangle(270)
xticks([]);
henry.TickDir='out';
set(henry,'FontName','Linux Libertine Display O','FontSize',22,'yscale','log');
title('Inheritance','FontSize',30,'FontName','Linux Libertine Display O');












% DATA NECESSARY FOR ADDITIONAL FIGURES
restot1=sum(bcs1(:,3));
restot2=sum(bcs2(:,3));
resfra1=bcs1(:,3)/restot1;
resfra2=bcs2(:,3)/restot2;
rescum1=zeros(size(resfra1));
rescum2=zeros(size(resfra2));
for cterant=2:length(resfra1)
    rescum1(cterant)=rescum1(cterant-1)+resfra1(cterant);
end
for dterant=2:length(resfra2)
    rescum2(dterant)=rescum2(dterant-1)+resfra1(dterant);
end










f1=figure; % COMPOSITITION FUNCTION
a1=axes('Position',[0.129568106312292 0.11 0.775431893687706 0.784806924101198]);
plot((1:length(rescum1))/length(resfra1),rescum1,'b','DisplayName','No inheritance')
hold on
plot((1:length(rescum2))/length(resfra2),rescum2,'r','DisplayName','Inheritance')
hold on
line([0 0.1],[rescum1(round(length(rescum1)/10)) rescum1(round(length(rescum1)/10))],'Color','blue','LineStyle','--');
line([0.1 0.1],[0 rescum1(round(length(rescum1)/10))],'Color','blue','LineStyle','--');
line([0 0.1],[rescum2(round(length(rescum2)/10)) rescum2(round(length(rescum2)/10))],'Color','red','LineStyle','--');
line([0.1 0.1],[0 rescum2(round(length(rescum2)/10))],'Color','red','LineStyle','--');
line([0 0.2],[rescum1(round(length(rescum1)/5)) rescum1(round(length(rescum1)/5))],'Color','blue','LineStyle','-.');
line([0.2 0.2],[0 rescum1(round(length(rescum1)/5))],'Color','blue','LineStyle','-.');
line([0 0.2],[rescum2(round(length(rescum2)/5)) rescum2(round(length(rescum2)/5))],'Color','red','LineStyle','-.');
line([0.2 0.2],[0 rescum2(round(length(rescum2)/5))],'Color','red','LineStyle','-.');
ylabel('Fraction of response (% of cells)','FontSize',22,'FontName','Linux Libertine Display O')
xlabel('Fraction of response (% of families)','FontSize',22,'FontName','Linux Libertine Display O')
set(a1,'FontName','Linux Libertine Display O','FontSize',22)
title('Response composition','FontSize',30,'FontName','Linux Libertine Display O');
legend1 = legend(a1,{'No Inheritance','Inheritance'});
set(legend1,'Position',[0.648947951273533 0.743009320905459 0.240310077519379 0.134487350199734]);
xticks([0 0.25 0.5 0.75 1])
xticklabels([0 25 50 75 100])
yticks([0 0.25 0.5 0.75 1])
yticklabels([0 25 50 75 100])










% EXPONENTIAL GROWTH PLOT
tots=zeros(6,2);
for yterant=1:4
    tots(yterant+1,1)=sum(bcd1(:,33+yterant));
    tots(yterant+1,2)=sum(bcd2(:,33+yterant));
end
tots(1,1)=length(bcd1(:,3));
tots(1,2)=length(bcd2(:,3));
tots(6,1)=sum(bcd1(:,3));
tots(6,2)=sum(bcd2(:,3));

famnum1=length(bcd1(:,3));
famnum2=length(bcd2(:,3));
cells1=zeros(5,50);
cells2=zeros(5,50);

for gterant=1:50
    rnum1=randi([1 famnum1],1,500);
    rnum2=randi([1 famnum2],1,500);
    for hterant=1:4
        cells1(hterant,gterant)=sum(bcd1(rnum1,33+hterant));
        cells2(hterant,gterant)=sum(bcd2(rnum2,33+hterant));
    end
    cells1(5,gterant)=sum(bcd1(rnum1,3));
    cells2(5,gterant)=sum(bcd2(rnum2,3));
end










figure; % GROWTH EXPONENTIAL PLOTS
a1=subplot(1,2,1)
boxplot(cells1','Labels',{'1','2','3','4','5'})
xlabel('Days')
ylabel('Cells')
title('No Inheritance')
hold on
plot(a1,0.5:0.1:5.5,500*exp(0.7864*[0.5:0.1:5.5]),'b-.')
%plot(a1,0.5:0.1:5.5,500*exp(0.7838*[0.5:0.1:5.5]),'b:')
%plot(a1,0.5:0.1:5.5,500*exp(0.7889*[0.5:0.1:5.5]),'b:')
set(a1,'FontName','Linux Libertine Display O','FontSize',22)
a2=subplot(1,2,2)
boxplot(cells2','Labels',{'1','2','3','4','5'})
xlabel('Days')
ylabel('Cells')
title('Inheritance')
hold on
plot(a2,0.5:0.1:5.5,500*exp(0.8345*[0.5:0.1:5.5]),'r-.')
%plot(a2,0.5:0.1:5.5,500*exp(0.8345*[0.5:0.1:5.5]),'r:')
%plot(a2,0.5:0.1:5.5,500*exp(0.8345*[0.5:0.1:5.5]),'r:')
set(a2,'FontName','Linux Libertine Display O','FontSize',22)













% CALCULATING GROWTH RATES
rates1=bsxfun(@minus,log(dbcd1(:,[34:37 3])),...
        [log(ones(size(dbcd1(:,1),1),1))...
        log(dbcd1(:,34:37))]);
rates2=bsxfun(@minus,log(dbcd2(:,[34:37 3])),...
        [log(ones(size(dbcd2(:,1),1),1))...
        log(dbcd2(:,34:37))]);

    
    
    

    
    
    
    

    
    
figure; % PLOT GROWTH RATE VARIANCES BETWEEN DAYS
a1=subplot(1,4,1)
boxplot(var(rates1'),'Labels',{'Time'})
ylim([0 2])
ylabel('Variance')
set(a1,'FontName','Linux Libertine Display O','FontSize',22)
a2=subplot(1,4,3)
boxplot(var(rates2'),'Labels',{'Time'})
ylim([0 2])
set(a2,'FontName','Linux Libertine Display O','FontSize',22)
a3=subplot(1,4,2)
boxplot(var(rates1),'Labels',{'Family'})
ylim([0 2])
set(a3,'FontName','Linux Libertine Display O','FontSize',22)
a4=subplot(1,4,4)
boxplot(var(rates2),'Labels',{'Family'})
ylim([0 2])
set(a4,'FontName','Linux Libertine Display O','FontSize',22)
mean(rates1(:))
mean(rates2(:))













namesvec2=[...
        "\color{gray} Simulation identifier",...
        "\color{gray} Barcode",...
        "\color{blue} Family size",...
        "\color{red} Competitors",...
        "\color{red} Naive cell's first division clock",...
        "\color{red} Naive cell's division clock",...
        "\color{red} Contact requirement to become activated",...
        "\color{red} Myc value upon becoming activated",...
        "\color{red} APC signalling strength",...
        "Time to naive cell activation",...
        "Contact received upon naive cell activation",...
        "Myc increase upon naive cell activation",...
        "Myc clock value upon naive cell activation",...
        "Time to first break between naive cell and APC",...
        "Contact received at first break between naive cell and APC",...
        "Myc increase at first break between naive cell and APC",...
        "Myc clock value at first break between naive cell and APC",...
        "Time to first division",...
        "Total contact received by naive cell",...
        "Total Myc increace received by naive cell",...
        "Myc clock value at first division",...
        "Time to first break between naive cell and APC after activation",...
        "Contact received upon first break between naive cell and APC after activation",...
        "Myc increase received upon first break between naive cell and APC after activation",...
        "Myc clock value upon first break between naive cell and APC after activation",...
        "Contact events in naive cell",...
        "\color{red} Volume per cell",...
        "\color{gray} Total number of cells in simulation space",...
        "\color{red} Number of CD8+ cells competing for APC contact",...
        "\color{red} Number of reticular cells in simulation space",...
        "Maximal contact received in family",...
        "Total contact events in family",...
        "Active daughter cells at simulation end",...
        "Family size at day 1.",...
        "Family size at day 2.",...
        "Family size at day 3.",...
        "Family size at day 4.",...
        "Family size at day 5."...
        ];

corrmat1=NaN(32,32);
corrmat2=NaN(32,32);

for qq=1:2

    if qq==1
        bcd=bcd1;
    else
        bcd=bcd2;
    end 
    for pp=1:32
        for nn=1:32

            list=min(isfinite(bcd(:,pp)),isfinite(bcd(:,nn)));
            variables=bcd(list,pp);
            measures=bcd(list,nn);
            xm=mean(variables(:));
            ym=mean(measures(:));
            corr=sum((variables(:)-xm).*(measures(:)-ym))...
                     /sqrt(sum((variables(:)-xm).^2)...
                     *sum((measures(:)-ym).^2));
            if qq==1
                corrmat1(pp,nn)=corr;
            else
                corrmat2(pp,nn)=corr;
            end
        end
    end
end

sum(isnan(corrmat1(:)));
sum(isnan(corrmat2(:)));

trans =[6 4 29 27 30 5 7 8 9 11 18 10 23 19 26 31 15 14 17 16 22 25 13 12 24 20 21 32 3];
num=length(trans);
corrmat3=NaN(num,num);

for ii = 1:num
    for jj = 1:num
        corrmat3(ii,jj)=(corrmat1(trans(jj),trans(ii)));
    end
end

namesvec3=namesvec2;
for zterant=1:length(trans)
    tran=trans(zterant);
    namesvec3(tran)=namesvec2(tran)+"  "+num2str(zterant);
end

figure;
colormap('redblue');
leonardo=subplot(3,1,1);
imagesc(corrmat3(1:9,10:num-1))
set(leonardo,'xdir','reverse');
yticks(1:9);
yticklabels(namesvec3(trans(1:9)));
xticks([])
title("Correlation between prediction factors");
caxis([-1 1]);
set(gca,'TickLength',[0 0])

donatello=subplot(3,1,2);
imagesc(corrmat3(10:num-1,10:num-1));
xticks([]);
yticks(1:num-10)
yticklabels(namesvec3(trans(10:num-1)))
set(donatello,'xdir','reverse')
colorbar;
caxis([-1 1]);
set(gca,'TickLength',[0 0])

michelangelo=subplot(3,1,3);
imagesc(corrmat3(num,10:num-1));
xticks(1:num-10);
xticklabels(10:num-1);
yticks(1);
yticklabels(namesvec3(trans(num)));
caxis([-1 1])
set(michelangelo,'xdir','reverse');
set(gca,'TickLength',[0 0])


