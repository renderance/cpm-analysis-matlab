function barcodemat = barcode_family_extractor(bcmat,folder,maxtime,graph)

    fprintf(1,'Working to extract interesting family properties.\n');

    factarray=[];

    for jj =1:length(bcmat(:,2))
       
        [maxcontact,contactevents,activetermini,day1,day2,day3,day4,day5]=descendanttracker2(bcmat(jj,2),bcmat(jj,2),graph,folder,maxtime,0);
        factarray=cat(1,factarray,[maxcontact,contactevents,activetermini,day1,day2,day3,day4,day5]);
        
    end
    
    barcodemat=cat(2,bcmat,factarray);

end

function [maxcontact_o, contactevents_o, activetermini_o, day_1_o, day_2_o, day_3_o, day_4_o, day_5_o] =...
    descendanttracker2(mothercode,barcode,graph,folder,maxtime,nest)

    fprintf(1,strcat('\tWorking on:\tbarcode:\t',num2str(mothercode),'\tdescendant:\t',num2str(barcode),'\tlayer:\t',num2str(nest),'\n'));

    filename=strcat(folder,'/cell_',num2str(barcode),'.csv');
    
    temp=importdata(filename)';
    
    timet_i=temp.data(end,1);
    times_i=temp.data(1,1);
    typet_i=temp.data(end,3);
    contact_list=temp.data(:,11);
    maxcontact_i=temp.data(end,10);

    % Count number of contact events in this cell.
    
    contactevents_i=0;
    prevstat=0;

    for jj=1:length(contact_list)

        newstat=contact_list(jj);

        if prevstat==0 && newstat==1

            contactevents_i=contactevents_i+1;

        end

        prevstat=newstat;

    end
    
    % Check if cell is terminal and active.
    
    if typet_i == 6 && timet_i==maxtime
        
        activetermini_i=1;
        
    else
        
        activetermini_i=0;
        
    end
    
    % Check if cell was present at end of day 1.
    
    if times_i<=14400 && timet_i>=14400
        
        day_1_i=1;
        
    else
        
        day_1_i=0;
        
    end
    
    % Check if cell was present at end of day 2.
    
    if times_i<=28800 && timet_i>=28800
        
        day_2_i=1;
        
    else
        
        day_2_i=0;
        
    end
    
    % Check if cell was present at end of day 3.
    
    if times_i<=43200 && timet_i>=43200
        
        day_3_i=1;
        
    else
        
        day_3_i=0;
        
    end
    
    % Check if cell was present at end of day 4.
    
    if times_i<=57600 && timet_i>=57600
        
        day_4_i=1;
        
    else
        
        day_4_i=0;
        
    end
    
    % Check if cell was present at end of day 5.
    
    if times_i<=72000 && timet_i>=72000
        
        day_5_i=1;
        
    else
        
        day_5_i=0;
        
    end
    
    % Reconcile with daughter data.
    
    for jj=1:length(graph.Nodes.Name)
            
        if compose("%s",graph.Nodes.Name{jj})==compose("%s",num2str(barcode))

            descendants = successors(graph,graph.Nodes.Name{jj});

        end

    end
    
    if ~isempty(descendants)

        for ii = 1:length(descendants)

            descendant=descendants{ii};
            [maxcontact_d, contactevents_d, activetermini_d, day_1_d, day_2_d, day_3_d, day_4_d, day_5_d]=...
                descendanttracker2(mothercode,descendant,graph,folder,maxtime,nest+1);
            
            maxcontact_i=max(maxcontact_d,maxcontact_i);
            contactevents_i=contactevents_d+contactevents_i;
            activetermini_i=activetermini_d+activetermini_i;
            day_1_i=day_1_i+day_1_d;
            day_2_i=day_2_i+day_2_d;
            day_3_i=day_3_i+day_3_d;
            day_4_i=day_4_i+day_4_d;
            day_5_i=day_5_i+day_5_d;
            
        end
        
    end
        
    maxcontact_o=maxcontact_i;
    contactevents_o=contactevents_i;
    activetermini_o=activetermini_i;
    day_1_o=day_1_i;
    day_2_o=day_2_i;
    day_3_o=day_3_i;
    day_4_o=day_4_i;
    day_5_o=day_5_i;

end