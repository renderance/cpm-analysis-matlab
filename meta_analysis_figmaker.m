
    fprintf(1,'Making great figures!\n');

    %% Import and sort all the data! 
    barcodedata=importdata('/home/glenn/Documents/morpheus/sim_outputs2/analysis_meta_barcode_data.csv');
    outputfolder='/home/glenn/Documents/morpheus/sim_outputs2/';

    barcodesorted=sortrows(barcodedata,[3 33],'descend');



    %% Make the plot that shows the variation in family size outcome!
    humphrey=figure('visible','off','Position', [10 10 1200 900]);
    timothy=subplot(1,2,1);
    bar(timothy,barcodesorted(:,3),'b');
    hold on
    bar(timothy,barcodesorted(:,33),'r');
    ylabel(timothy,'Family size / Number of active daughters')
    xlabel(timothy,'Barcodes')
    xtickangle(270)
    xticks([]);
    timothy.TickDir='out';
    saveas(humphrey,strcat(outputfolder,'analysis_meta_1_barcode_ranking.png'));
    close(humphrey)
    clear humphrey

    humphrey=figure('visible','off','Position', [10 10 1200 900]);
    timothy=axes;
    bar(timothy,barcodesorted(:,3),'b');
    hold on
    bar(timothy,barcodesorted(:,33),'r');
    ylabel(timothy,'Family size / Number of active daughters')
    xlabel(timothy,'Barcodes')
    xtickangle(270)
    xticks([]);
    set(timothy,'YScale','log')
    timothy.TickDir='out';
    saveas(humphrey,strcat(outputfolder,'analysis_meta_1_barcode_ranking_logscale.png'));
    close(humphrey)
    clear humphrey






    %% Plot time to first division!
    prunedbarcodes=barcodesorted(barcodesorted(:,3)>1,:);

    james=figure('visible','off','Position', [10 10 1200 900]);
    dylan=axes;
    john=histogram(dylan,prunedbarcodes(:,18),100);
    xticks(0:14400:72000);
    xticklabels(0:1:5);
    ylabel(dylan,'Frequency')
    xlabel(dylan,'Time to first division (days)')
    xlim([0 72000])
    saveas(james,strcat(outputfolder,'analysis_meta_2_time_to_first_division.png'));
    close(james)
    clear james





    %% Plot time to activation!
    wally=figure('visible','off','Position', [10 10 1200 900]);
    henry=axes;
    janet=histogram(henry,prunedbarcodes(:,10),100);
    xticks(0:14400:72000);
    xticklabels(0:1:5);
    ylabel(henry,'Frequency')
    xlabel(henry,'Time to activation (days)')
    xlim([0 72000])
    saveas(wally,strcat(outputfolder,'analysis_meta_2_time_to_activation.png'));
    close(wally)
    clear wally



    %% Plot time between activation and break from APC!
    jeff=figure('visible','off','Position', [10 10 1200 900]);
    amelia=axes;
    anette=histogram(amelia,prunedbarcodes(:,22),100);
    xticks(0:600:3600)
    xticklabels(0:1:6);
    xlim([0 3600])
    ylabel(amelia,'Frequency')
    xlabel(amelia,'Time between activation and break away from APC (hours)')
    saveas(jeff,strcat(outputfolder,'analysis_meta_2_time_of_additional_stimulation.png'));
    close(jeff)
    clear jeff





    %% Plot family expansion!
    molly=figure('visible','off','Position', [10 10 2000 900]);
    jeffrey=subplot(1,2,1);
    alison=plot(jeffrey,1:5,barcodedata(:,34:38),'-o');
    xticks(1:5)
    xlabel(jeffrey,'Days')
    ylabel(jeffrey,'Family size')
    milton=subplot(1,2,2);
    gerome=boxplot(milton,barcodedata(:,34:38),'Labels',{'1','2','3','4','5'});
    xlabel(milton,'Days')
    ylabel(milton,'Family size')
    saveas(molly,strcat(outputfolder,'analysis_meta_3_family_expansion.png'));
    close(molly)
    clear molly




    %% Calculate the growth rates per day!
    barcodedata(barcodedata(:,34)<1,34)=1;    
    rates1=bsxfun(@minus,log(barcodedata(:,34:38)),...
        [log(ones(size(barcodedata(:,1),1),1))...
        log(barcodedata(:,34:37))]);

    hulio=figure('visible','off','Position', [10 10 2000 900]);
    quintin=subplot(1,2,2);
    danika=boxplot(quintin,rates1(:,:),'Labels',{'0-1','1-2','2-3','3-4','4-5'});
    ylabel(quintin,'Growth rates')
    xlabel(quintin,'Days')
    arnold=subplot(1,2,1);
    sarah=plot(arnold,0:4,rates1(:,:),'-o');
    xticklabels({'0-1','1-2','2-3','3-4','4-5'});
    xticks(0:4)
    ylabel(arnold,'Growth rates')
    xlabel(arnold,'Days')
    saveas(hulio,strcat(outputfolder,'analysis_meta_4_daily_growth_rate_means.png'));
    close(hulio)
    clear hulio






    %% Make a bunch of vectors for names and stuff!    
    namesvec=[...
        "Simulation identifier",...
        "Barcode",...
        "Family size",...
        "Competitors",...
        "Naive cell's first division clock",...
        "Naive cell's division clock",...
        "Contact requirement to become activated",...
        "Myc value upon becoming activated",...
        "APC signalling strength",...
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
        "Volume per cell",...
        "Total number of cells in simulation space",...
        "Number of CD8+ cells competing for APC contact",...
        "Number of reticular cells in simulation space",...
        "Maximal contact received in family",...
        "Total contact events in family",...
        "Active daughter cells at simulation end",...
        "Family size at day 1.",...
        "Family size at day 2.",...
        "Family size at day 3.",...
        "Family size at day 4.",...
        "Family size at day 5."...
        ];

    labelvec=[...
        "Family size at day 5.",...
        "Mean of Growth Rate per day.",...
        "Mean of Growth Rate per two days.",...
        "Mean of Growth Rate per three days.",...
        "Mean of Growth Rate per four days.",...
        "Growth rate per five days."...
        ];

    titlevec=[...
        "family_size_",...
        "growth_rate_1_",...
        "growth_rate_2_",...
        "growth_rate_3_",...
        "growth_rate_4_",...
        "growth_rate_5_"...
        ];



    %% Prepare a bunch of arrays to take means of the growth rates!
    mrates1=NaN(size(rates1,1),1);

    for zz=1:size(rates1,1)

        mrates1(zz)=mean(rates1(zz,:));

        if isnan(mrates1(zz))

            fprintf(1,'NaN mean detected at zz=%d\n',zz)

        end

    end
     
     
     
     
     
     
     
     
     
    %% Plot Myc at first division!
    measures=barcodedata(~isnan(barcodedata(:,21)),3);
    variables=barcodedata(~isnan(barcodedata(:,21)),21);
    jj=1;
    ii=21;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    damian=scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_family_vs_myc_at_first_division','.png'));
    close(walter)
    clear walter

    measures=mrates1(min(~isnan(barcodedata(:,21)),isfinite(mrates1(:))));
    variables=barcodedata(min(~isnan(barcodedata(:,21)),isfinite(mrates1(:))),21);
    jj=2;
    ii=21;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    donald=scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_mrates1_vs_myc_at_first_division','.png'));
    close(walter)
    clear walter









    %% Plot the maximal contact time in a family
    measures=barcodedata(~isnan(barcodedata(:,31)),3);
    variables=barcodedata(~isnan(barcodedata(:,31)),31);
    jj=1;
    ii=31;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    lawrence=scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_family_vs_maximal_contact_time','.png'));
    close(walter)
    clear walter

    measures=mrates1(min(~isnan(barcodedata(:,31)),isfinite(mrates1(:))));
    variables=barcodedata(min(~isnan(barcodedata(:,31)),isfinite(mrates1(:))),31);
    jj=2;
    ii=31;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    jacky=scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_mrates1_vs_maximal_contact_time','.png'));
    close(walter)
    clear walter







    %% Plot contact events in family!
    measures=barcodedata(~isnan(barcodedata(:,32)),3);
    variables=barcodedata(~isnan(barcodedata(:,32)),32);
    jj=1;
    ii=32;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    mohammed=scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_family_vs_contact_events_in_family','.png'));
    close(walter)
    clear walter

    measures=mrates1(min(~isnan(barcodedata(:,32)),isfinite(mrates1(:))));
    variables=barcodedata(min(~isnan(barcodedata(:,32)),isfinite(mrates1(:))),32);
    jj=2;
    ii=32;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_mrates1_vs_contact_events_in_family','.png'));
    close(walter)
    clear walter






    %% Plot contact events in naive!
    measures=barcodedata(~isnan(barcodedata(:,26)),3);
    variables=barcodedata(~isnan(barcodedata(:,26)),26);
    jj=1;
    ii=26;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),cat(2,'Intercept: ',...
              num2str(fit(2))),cat(2,'Pearson: ',num2str(corr))},...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_family_vs_contact_events_in_naive','.png'));
    close(walter)
    clear walter

    measures=mrates1(min(~isnan(barcodedata(:,26)),isfinite(mrates1(:))));
    variables=barcodedata(min(~isnan(barcodedata(:,26)),isfinite(mrates1(:))),26);
    jj=2;
    ii=26;

    walter=figure('visible','off','Position', [10 10 1200 900]);
    thomas=axes;
    scatter(thomas,variables(:),measures(:));
    ylabel(labelvec(jj))
    xlabel(namesvec(ii))
        reg=lsline(thomas);
        reg.Color='r';
        fit=polyfit(variables(:),measures(:),1);
        xm=mean(variables(:));
        ym=mean(measures(:));
        corr=sum((variables(:)-xm).*(measures(:)-ym))...
              /sqrt(sum((variables(:)-xm).^2)...
              *sum((measures(:)-ym).^2));
        text(thomas,...
            (thomas.XLim(1)+(thomas.XLim(2)-thomas.XLim(1))/2),...
            (thomas.YLim(1)+(thomas.YLim(2)-thomas.YLim(1))*9/10),...
            {cat(2,'Slope: ',num2str(fit(1))),...
             cat(2,'Intercept: ',num2str(fit(2))),...
             cat(2,'Pearson: ',num2str(corr))...
            },...
            'Color','r'...
            )
    saveas(walter,strcat(outputfolder,'analysis_meta_5_mrates1_vs_contact_events_in_naive','.png'));
    close(walter)
    clear walter








    %% Calculate the fractions of the response are made up by several cells!
    response=sum(barcodesorted(:,3));
    response_frac=barcodesorted(:,3)/response;
    response_cumul=zeros(size(response_frac));
    response_cumul(1)=response_frac(1);

    for cc=2:length(response_frac)

        response_cumul(cc)=response_cumul(cc-1)+response_frac(cc);

    end

    
    
    
    
    
    %% Plot the accumulation of response!
    freddy=figure('visible','off','Position', [10 10 1200 900]);
    tamara=axes;
    plot(1:length(response_cumul),response_cumul)
    text(tamara,...
            (tamara.XLim(1)+(tamara.XLim(2)-tamara.XLim(1))/2),...
            (tamara.YLim(1)+(tamara.YLim(2)-tamara.YLim(1))*2/10),...
            {cat(2,'First Barcode: ',num2str(response_cumul(1))),...
             cat(2,'10% of Barcodes: ',num2str(response_cumul(round(length(response_cumul)/10)))),...
             cat(2,'20% of Barcodes: ',num2str(response_cumul(round(length(response_cumul)/5))))...
            },...
            'Color','k'...
            )
    hold on
    rectangle('Position',[0 0 round(length(response_cumul)/10) response_cumul(round(length(response_cumul)/10))])
    rectangle('Position',[0 0 round(length(response_cumul)/5) response_cumul(round(length(response_cumul)/5))])
    saveas(freddy,strcat(outputfolder,'analysis_meta_5_fraction_of_response','.png'));
    close(freddy) 
    clear freddy










    %% Write the correlations between interesting factors to a file!
    fid=fopen(strcat(outputfolder,'analysis_meta_5_factor_correlations.csv'),'w');

    for pp=1:32

        for nn=1:32


            list=min(isfinite(barcodedata(:,pp)),isfinite(barcodedata(:,nn)));
            variables=barcodedata(list,pp);
            measures=barcodedata(list,nn);
            xm=mean(variables(:));
            ym=mean(measures(:));
            corr=sum((variables(:)-xm).*(measures(:)-ym))...
                     /sqrt(sum((variables(:)-xm).^2)...
                     *sum((measures(:)-ym).^2));
            fprintf(fid,strcat('%d\t',num2str(pp),'\t',num2str(nn),'\n'),corr);


        end

    end

    fclose(fid);
    corrlist=readtable(strcat(outputfolder,'analysis_meta_5_factor_correlations.csv'));
    corrmat=reshape(corrlist.Var1,[32 32]);
    
    
    
    
    
    
    
    
    
%     %% Colormap that shit!
%     colormap('jet');
%     gwello=figure('visible','off','Position', [10 10 1200 1000]);
%     imagesc(corrmat);
%     xticks(1:31)
%     xtickangle(270)
%     xticklabels(namesvec)
%     yticks(1:31)
%     yticklabels(namesvec)
%     colorbar;
%     title("Correlation between prediction factors");
%     saveas(gwello,strcat(outputfolder,'analysis_meta_5_factor_correlations','.png'));
%     close(gwello)
%     clear gwello








    %% Update the namesvec with colors!
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

    
    
    
    
    %% Prepare and remake the correlation matrix!
    corrmat2=zeros(27);
    
    %trans=[5 6 4 27 14 15 16 17 10 18 7 11 19 23 8 9 12 13 24 20 25 22 26 21 31 32 3];
    trans =[5 6 4 27 8 9 7 11 23 19 31 22 25 20 12 24 13 16 17 14 15 18 10 26 21 32 3];
    
    for ii = 1:27

        for jj = 1:27

            corrmat2(ii,jj)=abs(corrmat(trans(jj),trans(ii)));

        end

    end

    
    
    %% Plot the correlation matrix nicely!
    
    gwello=figure('visible','off','Position', [10 10 1200 1000]);
    colormap('hot');
    donatello=axes;
    imagesc(corrmat2);
    xticks(1:27)
    xtickangle(270)
    xticklabels(namesvec2(trans))
    yticks(1:27)
    yticklabels(namesvec2(trans))
    set(donatello,'xdir','reverse')
    colorbar;
    caxis([0 1]);
    title("Correlation between prediction factors");
    saveas(gwello,strcat(outputfolder,'analysis_meta_5_factor_correlations_sorted_2','.png'));
    close(gwello)
    clear gwello