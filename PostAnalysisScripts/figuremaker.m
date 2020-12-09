
%% Make image of all found events (white) and grouped locations (red)
imagesc(sum(Data1,3))
axis image
hold on
for i=1:size(LocatStore,2);
    plot(LocatStore(1,i).PSFfinal(:,2),LocatStore(1,i).PSFfinal(:,1),'wx')
end

for i=1:size(GroupLocat,2);
    plot(GroupLocat(1,i).Centroid(1,2),GroupLocat(1,i).Centroid(1,1),'ro')
end

%% Make image of ensemble cumulative desorption and adsorption times
addpath('U:\Lydia\MATLAB\CumulDist\')

[dd,id]=cumuldist(Ensemble.Dwell,unique(Ensemble.Dwell));

% % fit the data
% model_fun = 'a*exp(-x/t)+46.5'; %change 46.5 accordingly (e.dataSpace/2)+e.dataTime
% coefs = {'a','t'};
% exp_model = fittype(model_fun,'coefficients',coefs,'independent','x');
% spt=46.5*5;%change 46.5 accordingly (e.dataSpace/2)+e.dataTime
% fitop=fitoptions('Lower',[0 0],'Upper',[2,inf],'StartPoint',[1 spt]);
% [cfun,gof] = fit(dd,id,exp_model,'options',fitop);  % no quotes for a fittype model!
% model_coefs = coeffvalues(cfun);
% model_fit = model_coefs(1)*exp(model_coefs(2)*dd)+46.5;%change 46.5 accordingly (e.dataSpace/2)+e.dataTime

% plot the data
figure
plot(dd,id)
set(gca,'YScale','log','FontSize',14);
title('Desorption ensemble cumulative distribution');
ylabel('P(t">t)');xlabel('Time (ms)')
hold on
% plot(dd,model_fit,'--')
% text(strcat('td = ',num2str(model_coefs(2)),' ms; R^2 = ',num2str(gof.rsquare)))

% do for association data
[da,ia]=cumuldist(Ensemble.Assoc,unique(Ensemble.Assoc));
figure
plot(da,ia)
set(gca,'YScale','log','FontSize',14);
title('Association ensemble cumulative distribution');

%% Make super-resolution figure
ylabel('P(t">t)');xlabel('Time (ms)')

%plot all individual site assoc distributions
figure
hold on
set(gca,'YScale','log','FontSize',14);
for i=1:size(GroupLocat,2);
	if isempty(GroupLocat(1,i).Assoc)==0
        [dd,id]=cumuldist(GroupLocat(1,i).Assoc,unique(GroupLocat(1,i).Assoc));    
        plot(dd,id)
    end
end
