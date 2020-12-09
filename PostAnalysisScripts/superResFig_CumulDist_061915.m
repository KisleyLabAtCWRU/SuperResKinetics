% e.nzoom=6;
close all
clear all


filelist=1:32;
for i=1:numel(filelist)
    
    %load the data
    fnum=num2str(filelist(i));
    load(strcat(fnum,'_analyzed.mat'))
    
    %create super res image
    e.nzoom=6;
    addpath('U:\Lydia\MATLAB\SuperPosition-20130215\LRG_SuperRes_Kinetics_Final')
    superdata=LRG_SuperRes_GenerateSR(LocatStore,e);
    addpath('U:\Lydia\MATLAB\CumulDist\')

    [dd,id]=cumuldist(Ensemble.Dwell,unique(Ensemble.Dwell));



    % plot the data
    figure
    subplot(1,2,1)
    plot(dd,id,'ko')
    set(gca,'YScale','log','FontSize',14);
    title('Desorption ensemble cumulative distribution');
    ylabel('P(t">t)');xlabel('Time (ms)')
    hold on
    axis square

    subplot(1,2,2)
    imagesc(superdata)
    axis image
    xlim([1000 2000])
    ylim([1000 2000])
    caxis([5 20])
    title(fnum)
    set(gca,'xtick',[],'ytick',[])

    savename=strcat(fnum,'_super_dwell');
    save(savename,'superdata','dd','id')
    
    clearvars -except filelist i
end
% 
% for i=1:size(Data1,3)
%     I(i)=sum(sum(Data1(:,:,i),1),2);
% end
% % plot(I)
% 
% movies=[];
% figure
% for i=100:600
%     imagesc(Data1(250:450,100:300,i))
%     colormap(gray)
%     caxis([0 1000])
%     axis image
%     set(gca,'xtick',[],'ytick',[])
%     movies=[movies getframe]; 
% end

