% Single site kinetic plots
% LK021814

% % Look at sites w/ more than n events
% indexgl=[]
% n=150;
% for i=1:size(GroupLocat,2)
%     if size(GroupLocat(1,i).Dwell,2)>n
%         indexgl=[indexgl,i];
%     end
% end

indexgl=[738]; %indeces from GroupLocat to plot

addpath('U:\Lydia\MATLAB\CumulDist\')

for i=1:numel(indexgl)

    j=indexgl(i);
    location=GroupLocat(1,j).Centroid;
    locattxt=num2str(location(1,1:2));
    
    [dd,id]=cumuldist(GroupLocat(1,j).Dwell,unique(GroupLocat(1,j).Dwell));
    [da,ia]=cumuldist(GroupLocat(1,j).Assoc,unique(GroupLocat(1,j).Assoc));
    
    figure
    subplot(1,2,1)
    plot(dd,id,'-ko','LineWidth',2,'MarkerSize',6)
    set(gca,'YScale','log','FontSize',14);
    title(strcat('Cumulative distribution for single site ',num2str(j),'-',locattxt));
    ylabel('Desorption P(t">t)','FontSize',14);xlabel('Time (ms)','FontSize',14)
    hold on

    subplot(1,2,2)
    plot(da,ia,'-ko','LineWidth',2,'MarkerSize',6)
    set(gca,'YScale','log','FontSize',14);
    ylabel('Association P(t">t)','FontSize',14);xlabel('Time (ms)','FontSize',14)
    hold on
end