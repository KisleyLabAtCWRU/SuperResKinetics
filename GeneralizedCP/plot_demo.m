function plot_demo(eff,s,eff_fit,MDL)
eff = eff(480:end);
s = s(3:6,480:end);
eff_fit = eff_fit(:,480:end);
l = numel(s(:,1));
figure
for i = 1:l
    subplot(l,1,i)
    plot(eff,'b:')
    hold on
    plot(s(i,:),'r','LineWidth',1.5)
    set(gca, 'xtick',[],'ytick',[])
    ylim([0.6,0.8])
    xlim([0,length(s)+1])
    hold off
end
figure
a = numel(eff_fit(:,1));
%indx = round(logspace(0,log10(a),6));
indx = [2,3,12];
for i = 1:3
    subplot(4,1,i)
    plot(eff,'b:')
    hold on
    plot(eff_fit(indx(4-i),:),'k','LineWidth',1.5)
    set(gca, 'xtick',[],'ytick',[])
    ylim([0.6,0.8])
    xlim([0,length(s)+1])
    hold off
end
i = 4;
subplot(4,1,i)
plot(MDL,'s-k')
set(gca,'ytick',[])