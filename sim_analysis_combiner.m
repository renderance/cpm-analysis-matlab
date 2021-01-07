function [xmean,ymean,xvar,yvar,covar,finstat,undiv,div,all,bcmat]=...
    sim_analysis_combiner(outputfolder,folder,cellfiles,timefile,identifier,targettime)
    % This function calls all analysis functions for one simulation.

    tic;
    
    % Make abundance plots of all populations for this simulation.
    %maxtime=cellcount_all_timeplot(outputfolder,folder,timefile,identifier);
    %cellcount_apc_timeplot(outputfolder,folder,timefile,identifier)
    %cellcount_cd8_timeplot(outputfolder,folder,timefile,identifier)
    %cellcount_cd4_timeplot(outputfolder,folder,timefile,identifier)
    
    % Check if simulation finished.
    %finstat=(maxtime~=targettime);
    
    % Load motility data.
    maxtime=72000;
    [xcos,ycos,famv,idsv,momv,~,~,~,~,motv]=motility_data_loader(cellfiles,maxtime);
    
    % Perform a full motility analysis for this simulation.
    %motility_all_traceplot(outputfolder,xcos,ycos,motv,identifier)
    %motility_subset_traceplot(outputfolder,xcos,ycos,cellfiles,identifier)
    %motility_alt_measure_displacement(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)
    %motility_mean_square_displacement(outputfolder,xcos,ycos,motv,cellfiles,maxtime,identifier)
    %motility_angle_distribution(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)
    %motility_turning_angle_distribution(outputfolder,xcos,ycos,cellfiles,maxtime,identifier)
    %[xmean,ymean,xvar,yvar,covar]=motility_bivariate_normal_fitter(outputfolder,xcos,ycos,maxtime,identifier);
    
    % Make a decendancy graph.
    [undiv,div,all,bcmat,graphobj]=decendancy_branch_grapher(outputfolder,folder,identifier,idsv,famv,momv);
    
    % Extract information to forward to meta-analysis.
    bcmat=barcode_naive_extractor(bcmat,folder,maxtime);
    bcmat=barcode_time_extractor(bcmat,folder,timefile);
    bcmat=barcode_family_extractor(bcmat,folder,maxtime,graphobj);
    
    timer=toc;
    fprintf(1,'Full analysis of simulation %d took %d.\n',identifier,timer);
    
end

