function [undivided,divided,total,barcodes,graph]=decendancy_branch_grapher(outputfolder,folder,identifier,idsv,famv,momv)

    fprintf(1,'Making decendancy graph.\n')

    graph=digraph;
    
    decendancy=[momv(famv~=0) idsv(famv~=0)];
    
    for ii=1:length(decendancy(:,1))
        
        graph=addedge(graph,num2str(decendancy(ii,1)),num2str(decendancy(ii,2)));
        
    end

    graph=rmnode(graph,num2str(0));
    
    humbolt=figure('visible','off','Position', [10 10 1600 1200]);
    title(['Family Tree of simulation ',identifier])
    corgi=plot(graph,'Layout','layered','direction','right');
    box(corgi.Parent,'off');
    corgi.Parent.XTick=[];
    corgi.Parent.YTick=[];
    labelnode(corgi,graph.Nodes{:,'Name'},graph.Nodes{:,'Name'});
    saveas(humbolt,strcat(outputfolder,'decendancy_graph_',num2str(identifier),'.png'));
    close(humbolt)
    
    undivided=0;
    divided=0;
    ancestorvec=strings;

    for jj=1:length(graph.Nodes.Name)

        if isempty(predecessors(graph,graph.Nodes.Name(jj)))

            if isempty(successors(graph,graph.Nodes.Name(jj)))

                undivided=undivided+1;

            else

                divided=divided+1;

            end

            ancestorvec=cat(1,ancestorvec,compose("%s",graph.Nodes.Name{jj}));

        end

    end

    ancestorvec=ancestorvec(2:end);
    total=undivided+divided;
    ancestors=strings(1,length(graph.Nodes.Name));

    for kk=1:length(graph.Nodes.Name)

        ancestors(kk)=findancestor(graph,graph.Nodes.Name(kk));

    end

    leafcount=zeros(size(ancestorvec));

    for ll=1:length(ancestors)

        if isempty(successors(graph,graph.Nodes.Name(ll)))

            leafcount(ancestorvec==ancestors(ll))=leafcount(ancestorvec==ancestors(ll))+1;

        end

    end

    barcodes=cat(2,repmat(identifier,length(ancestorvec),1),double(ancestorvec),leafcount);

end

function ancestor = findancestor(object,nodename)
        
    if isempty(predecessors(object,object.Nodes.Name(object.Nodes.Name==compose("%s",nodename{1}))))

        ancestor = nodename;

    else

        ancestor = findancestor(object,predecessors(object,...
            object.Nodes.Name(object.Nodes.Name==compose("%s",nodename{1}))...
        ));

    end

end