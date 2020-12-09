%% avgcumuldist
% LK 010714

%% Average cumulative distribution function
% load cumulative distributions from multiple files and calculate the
% average and standard deviation for a plot
clear all
close all

% User input
filelist=[27:31];
dwell=0;
makefigure=1;
savedata=1;

% Load data and calculate individual cumulative distributions
addpath('U:\Lydia\MATLAB\CumulDist\')
alltime=[];
for i=1:numel(filelist)
    fname=strcat(num2str(filelist(i)),'_analyzed.mat');
    load(fname,'Ensemble')
    if dwell==1
        [time,prob]=cumuldist(Ensemble.Dwell,unique(Ensemble.Dwell));
    else
        [time,prob]=cumuldist(Ensemble.Assoc,unique(Ensemble.Assoc));
    end
    alltime=[alltime;time];
    store(1,i).data=[time,prob];    
end

% Organize time values calculated - will be final x values
uniquetime=unique(alltime);
sorttime=sort(uniquetime);

% Calculate average and std
for i=1:numel(sorttime)
    probtoavg=[];
    for j=1:size(store,2)
        for k=1:size(store(1,j).data,1)
            if store(1,j).data(k,1)==sorttime(i)
                probtoavg=[probtoavg,store(1,j).data(k,2)];
            end
        end    
    end
    meanprob=mean(probtoavg);
    stdprob=std(probtoavg);   
    final(i,:)=[sorttime(i),meanprob,stdprob];
end

% make figure
if makefigure==1
    figure
    errorbar(final(:,1),final(:,2),final(:,3),'ko','LineWidth',2);
    set(gca,'FontSize',14,'yscale','log')
    xlabel('Time (ms)','FontSize',14)
    ylabel('P(t`>t)','FontSize',14)
end

% save data
if savedata==1
    if dwell==1
        d='dwell';
    else
        d='assoc';
    end
    savename=strcat(num2str(filelist(1)),'_',num2str(filelist(numel(filelist))),'avgcumul_',d,'.mat');
    save(savename,'final')
end

